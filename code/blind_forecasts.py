"""Frozen blind-distance forecasts F1-F5 + local-exponent drift gate (idea L, SL3 pick 1).
Reads ONLY stored d=3,5,7 files (truth_d{3,5,7}_*, truthextra_d{3,5,7}_*, direct_d{3,5,7}_*). Never opens any d>=9 file.
numpy only. Usage: python3 blind_forecasts.py  -> results/BLIND_D_FORECASTS.json, results/BLIND_D_FORECASTS_gate.txt"""
import json, glob, os, sys, hashlib
import numpy as np
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','results')
DS=(3,5,7); PTAB=(1e-3,1e-4,1e-5,1e-6,2e-3,4e-3); P0=1e-3; FD=(9,11,13); NSIM=400000; SEED=20261125
def pstr(p): return f'{p:g}'
def load_ref():
    """reference per (d,p): IVW of IS truth (run1+run2 IVW) and direct k/n where a direct file exists. block rate, rounds=d."""
    T={}
    for f in glob.glob(f'{R}/truth_d*_p*.json'):
        x=json.load(open(f))
        if x['d'] in DS: T[(x['d'],x['p'])]=(x['truth'],x['se'])
    for f in glob.glob(f'{R}/truthextra_d*_p*.json'):
        y=json.load(open(f))
        if y['d'] not in DS: continue
        t,s=T[(y['d'],y['p'])]; w1,w2=1/s**2,1/y['se']**2
        T[(y['d'],y['p'])]=((w1*t+w2*y['truth'])/(w1+w2),(w1+w2)**-.5)
    ref={}
    for k,(t,s) in T.items():
        fn=f'{R}/direct_d{k[0]}_p{pstr(k[1])}.json'
        if os.path.exists(fn):
            x=json.load(open(fn)); dv=x['k']/x['n']; ds=x['k']**.5/x['n']
            w1,w2=1/s**2,1/ds**2; t,s=(w1*t+w2*dv)/(w1+w2),(w1+w2)**-.5
        ref[k]=(t,s)
    return ref
def eps_of(pl,d): return -np.expm1(np.log1p(-2*pl)/d)/2
def block_of(eps,d): return -np.expm1(d*np.log1p(-2*eps))/2
def ln_eps(ref,d,p):
    t,s=ref[(d,p)]; e=eps_of(t,d); de=(1-2*t)**(1/d-1)/d   # d eps / d pL
    return np.log(e), s*de/e
def mvn(mu,cov,rng,n): return rng.multivariate_normal(mu,cov+1e-300*np.eye(len(mu)),size=n,method='cholesky')
def summ(blk,eps_ln,hs):
    """blk, eps_ln: dict d-> arrays of draws"""
    out={}
    for d in FD:
        b=blk[d]; e=np.exp(eps_ln[d])
        out[f'd{d}']={'block_rate_median':float(np.median(b)),'block_rate_lo95':float(np.percentile(b,2.5)),'block_rate_hi95':float(np.percentile(b,97.5)),
          'per_round_median':float(np.median(e)),'per_round_lo95':float(np.percentile(e,2.5)),'per_round_hi95':float(np.percentile(e,97.5))}
    for d in FD:
        a=eps_ln[d-2] if d-2 in eps_ln else None
        if a is None: continue
        lam=np.exp(a-eps_ln[d])
        out[f'd{d}']['Lambda_from_prev']={'pair':f'{d-2}->{d}','median':float(np.median(lam)),'lo95':float(np.percentile(lam,2.5)),'hi95':float(np.percentile(lam,97.5))}
    return out
def linfit(X,y,sig):
    W=1/sig**2; A=X.T@(X*W[:,None]); cov=np.linalg.inv(A); beta=cov@(X.T@(W*y)); return beta,cov
def run_F1_F2_F3(ref,rng):
    ds=np.array(DS,float); y=[];s=[]
    for d in DS: a,b=ln_eps(ref,d,P0); y.append(a); s.append(b)
    y=np.array(y); s=np.array(s)
    X=np.c_[np.ones(3),ds]; beta,cov=linfit(X,y,s)
    def draws(beta,cov,extra=0.0):
        B=mvn(beta,cov,rng,NSIM); eps_ln={d:B[:,0]+B[:,1]*d for d in (7,)+FD}
        if extra>0: eps_ln={d:v+rng.normal(0,np.sqrt(extra),NSIM) if d in FD else v for d,v in eps_ln.items()}
        blk={d:block_of(np.exp(eps_ln[d]),d) for d in FD}; return blk,eps_ln
    res={}
    blk,el=draws(beta,cov); res['F1']=summ(blk,el,None)
    res['F1']['fit']={'ln_eps_intercept':float(beta[0]),'ln_eps_slope_per_d':float(beta[1]),'Lambda_per_distance_step_2':float(np.exp(-2*beta[1])),'chi2_dof1':float(((y-X@beta)**2/s**2).sum())}
    # F2: Lambda from d=5 to 7 only, geometric from d=7 anchor
    L=y[1]-y[2]; sL=np.sqrt(s[1]**2+s[2]**2)
    z5=rng.normal(y[1],s[1],NSIM); z7=rng.normal(y[2],s[2],NSIM); lam=z5-z7
    el={7:z7}
    for d in FD: el[d]=z7-lam*(d-7)/2
    blk={d:block_of(np.exp(el[d]),d) for d in FD}; res['F2']=summ(blk,el,None)
    res['F2']['fit']={'ln_Lambda_5to7':float(L),'se_ln_Lambda':float(sL)}
    # F3: F1 fit to d=3,5 only gives LOO residual at d=7; add square to the F1 variance at each forecast d
    slope=(y[1]-y[0])/2; pred7=y[1]+slope*2; v7=s[0]**2*0+0  # exact two-point line
    loo=float(y[2]-pred7)
    blk,el=draws(beta,cov,extra=loo**2); res['F3']=summ(blk,el,None); res['F3']['fit']={'loo_log_residual_d7':loo,'added_variance':loo**2}
    return res
def run_F4_F5(ref,rng):
    pts=[(d,p) for d in DS for p in PTAB]
    for k in pts: assert k in ref,k
    d=np.array([k[0] for k in pts],float); p=np.array([k[1] for k in pts]); h=(d+1)/2; x=p/1e-2
    yy=[];ss=[]
    for k in pts: a,b=ln_eps(ref,*k); yy.append(a); ss.append(b)
    yy=np.array(yy); ss=np.array(ss)
    # F4: ln eps - h ln p = sum_j (a_j + b_j h) x^j, j=0..2  (R4' form, coefficients linear in h)
    z=yy-h*np.log(p)
    def X4(h,x): return np.c_[np.ones_like(x),x,x**2,h,h*x,h*x**2]
    beta,cov=linfit(X4(h,x),z,ss)
    chi2=float((((z-X4(h,x)@beta)/ss)**2).sum()); dof=len(z)-6
    B=mvn(beta,cov,rng,NSIM); x0=P0/1e-2
    el={}
    for dd in (3,5,7)+FD:
        hh=(dd+1)/2; el[dd]=hh*np.log(P0)+B@X4(np.array([hh]),np.array([x0]))[0]
    blk={dd:block_of(np.exp(el[dd]),dd) for dd in FD}; r4=summ(blk,el,None); r4['fit']={'chi2':chi2,'dof':dof,'beta':beta.tolist()}
    # F5: R2-lin form  eps = p^h (c0(h) + c1(h) x),  c0=exp(a0+b0 h), c1=exp(a1+b1 h)*c0   (so c0+c1 x = c0 (1 + exp(a1+b1 h) x)); LM, numeric Jacobian
    lp=np.log(p)
    def f(th): return h*lp+th[0]+th[1]*h+np.log(np.maximum(1+(th[2]+th[3]*h)*x,1e-300))
    th=np.array([np.polyfit(h[p==P0],yy[p==P0]-h[p==P0]*np.log(P0),1)[1],np.polyfit(h[p==P0],yy[p==P0]-h[p==P0]*np.log(P0),1)[0],0.0,0.0])
    lam=1e-3
    def jac(th):
        J=np.empty((len(yy),4))
        for i in range(4):
            e=np.zeros(4); e[i]=1e-6; J[:,i]=(f(th+e)-f(th-e))/2e-6
        return J
    def chi(th):
        v=(((yy-f(th))/ss)**2).sum(); return v if np.isfinite(v) else np.inf
    cur=chi(th)
    for it in range(500):
        J=jac(th)/ss[:,None]; r=(yy-f(th))/ss; g=J.T@r; H=J.T@J
        step=np.linalg.lstsq(H+lam*np.diag(np.diag(H)+1e-12),g,rcond=1e-12)[0]; new=chi(th+step)
        if new<cur: th=th+step; lam=max(lam/3,1e-12); dec=cur-new; cur=new
        else:
            lam*=5; dec=1
            if lam>1e12: break
        if dec<1e-13 and it>20: break
    J=jac(th)/ss[:,None]; covth=np.linalg.pinv(J.T@J,rcond=1e-12)
    Bt=mvn(th,covth,rng,NSIM)
    el5={}
    for dd in (3,5,7)+FD:
        hh=(dd+1)/2; el5[dd]=hh*np.log(P0)+Bt[:,0]+Bt[:,1]*hh+np.log(np.maximum(1+(Bt[:,2]+Bt[:,3]*hh)*x0,1e-300))
    blk5={dd:block_of(np.exp(el5[dd]),dd) for dd in FD}; r5=summ(blk5,el5,None); r5['fit']={'chi2':float(cur),'dof':len(yy)-4,'theta':th.tolist(),'converged_iters':it}
    return r4,r5
def gate(ref,rng):
    lines=[]; cells=[]
    for p in (1e-3,1e-4,1e-5,1e-6):
        e=[ln_eps(ref,d,p) for d in DS]; y=np.array([a for a,b in e]); s=np.array([b for a,b in e])
        c=np.array([-1.,2.,-1.])   # ln drift = ln[Lambda(5->7)/Lambda(3->5)] = 2 ln e5 - ln e3 - ln e7
        D=float(c@y); sD=float(np.sqrt((c**2*s**2).sum())); T=D/sD
        # exact check: sigma=0 fixed-Lambda synthetic data gives drift 0 within 1e-12
        Lm=np.log(1.9); a0=np.log(4e-4); yf=np.array([a0-Lm*(d+1)/2 for d in DS]); exact=abs(float(c@yf))
        # calibrate: fixed-Lambda null (drift 0) with this cell's se, and 30% Lambda change alternative (ln 1.3 shift)
        null=rng.normal(0,1,(NSIM,3))*s; fa=float((np.abs((null@c)/sD)>1.96).mean())
        alt=rng.normal(0,1,(NSIM,3))*s; shift=np.log(1.3)
        pw=float(((np.abs(((alt@c)+shift)/sD)>1.96)|False).mean()); pw_dn=float((np.abs(((alt@c)-shift)/sD)>1.96).mean())
        usable=bool(fa<=0.075 and min(pw,pw_dn)>=0.80)
        cells.append({'p':p,'ln_drift':D,'drift_ratio_Lambda57_over_Lambda35':float(np.exp(D)),'se_ln_drift':sD,'T':T,'alarm_at_nominal_5pct':bool(abs(T)>1.96),
          'rel_se_refs':[float(v) for v in s],'false_alarm_rate':fa,'power_Lambda_up30':pw,'power_Lambda_down30':pw_dn,'exact_check_abs':exact,'exact_ok':bool(exact<1e-12),'usable':usable})
    return cells
if __name__=='__main__':
    rng=np.random.default_rng(SEED); ref=load_ref()
    out={'registered':'SL3 pick 1; see PLAN.md amendment BLIND-D FORECASTS','seed':SEED,'nsim':NSIM,'p':P0,'forecast_distances':list(FD),
         'inputs':{f'd{d}':{'reference_block_rate':ref[(d,P0)][0],'se':ref[(d,P0)][1]} for d in DS},
         'rounds_rule':'rounds=d; per-round eps=(1-(1-2pL)^(1/d))/2; block rate = (1-(1-2eps)^d)/2'}
    F=run_F1_F2_F3(ref,rng); F['F4'],F['F5']=run_F4_F5(ref,rng); out['forecasts']=F; out['drift_gate']=gate(ref,rng)
    out['code_sha256']=hashlib.sha256(open(os.path.abspath(__file__),'rb').read()).hexdigest()
    json.dump(out,open(f'{R}/BLIND_D_FORECASTS.json','w'),indent=1,sort_keys=True)
    print('wrote BLIND_D_FORECASTS.json')
