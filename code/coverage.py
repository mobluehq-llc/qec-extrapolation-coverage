# fits (R1..R5) on replicate datasets, delta-method 95% intervals, coverage vs truth, Clopper-Pearson.
import json,sys,numpy as np
from scipy.stats import beta
PS=np.array([4e-3,6e-3,8e-3,1e-2,1.5e-2]); X0=np.log(PS)
TARGETS=[1e-3,1e-4,1e-5,1e-6]
def wls(X,y,w):
    A=X.T@(X*w[:,None]); C=np.linalg.inv(A); b=C@(X.T@(w*y)); r=y-X@b; chi2=float((w*r*r).sum()); return b,C,chi2,len(y)-X.shape[1]
def fit_all(data,tp,ds=(3,5,7),scaled=False):
    """data[d]=(k[5],n[5]); returns {(recipe,d):(mu,var)} for log P at target p=tp"""
    xs=np.log(tp); out={}; Y={};Wt={}
    for d in data:
        k,n=data[d]; P=k/n; Y[d]=np.log(P); Wt[d]=k/(1-P)   # 1/var(log P)=k/(1-P)
    def sc(C,chi2,dof): return C*max(1.0,chi2/dof) if scaled and dof>0 else C
    for d in data:
        y,w=Y[d],Wt[d]
        X=np.c_[np.ones(5),X0]; b,C,c2,dof=wls(X,y,w); C=sc(C,c2,dof); g=np.array([1,xs]); out[('R1',d)]=(g@b,g@C@g)
        h=(d+1)/2; X=np.ones((5,1)); b,C,c2,dof=wls(X,y-h*X0,w); C=sc(C,c2,dof); out[('R2',d)]=(b[0]+h*xs,C[0,0])
        X=np.c_[np.ones(5),PS/1e-2,(PS/1e-2)**2]; b,C,c2,dof=wls(X,y-(d/2)*X0,w); C=sc(C,c2,dof); g=np.array([1,tp/1e-2,(tp/1e-2)**2]); out[('R4',d)]=(g@b+(d/2)*xs,g@C@g)
    # R3 shared
    Xs=[];ys=[];ws=[]
    for d in data:
        h=(d+1)/2; Xs.append(np.c_[np.ones(5),np.full(5,h)]); ys.append(Y[d]-h*X0); ws.append(Wt[d])
    X=np.vstack(Xs);y=np.concatenate(ys);w=np.concatenate(ws); b,C,c2,dof=wls(X,y,w); C=sc(C,c2,dof)
    for d in data:
        h=(d+1)/2; g=np.array([1,h]); out[('R3',d)]=(g@b+h*xs,g@C@g)
    # R5: Lambda from two smallest d (R1 per-d fits), extrapolate to larger d
    if 3 in data and 5 in data:
        m3,v3=out[('R1',3)]; m5,v5=out[('R1',5)]
        for d in data:
            if d>=7:
                s=(d-5)//2; out[('R5',d)]=(m5+s*(m5-m3),(1+s)**2*v5+s**2*v3)
    return out
def cp(k,n):
    lo=beta.ppf(.025,k,n-k+1) if k>0 else 0.0; hi=beta.ppf(.975,k+1,n-k) if k<n else 1.0; return lo,hi
def coverage(reps,truth,ds=(3,5,7),targets=TARGETS,scaled=False):
    """reps: list of dict d->(k,n) ; truth[(d,tp)]=p_L ; returns rows"""
    cnt={};tot={}
    for data in reps:
        for tp in targets:
            o=fit_all({d:data[d] for d in ds},tp,scaled=scaled)
            for (rc,d),(mu,var) in o.items():
                t=truth.get((d,tp))
                if t is None: continue
                hit=abs(np.log(t)-mu)<=1.96*np.sqrt(var)
                key=(rc,d,tp); cnt[key]=cnt.get(key,0)+int(hit); tot[key]=tot.get(key,0)+1
    rows=[]
    for key in sorted(cnt):
        k,n=cnt[key],tot[key]; lo,hi=cp(k,n); rows.append(dict(recipe=key[0],d=key[1],target=key[2],hits=k,reps=n,cov=k/n,cp_lo=lo,cp_hi=hi,
            miscal=int(k/n<0.90),miscal_firm=int(hi<0.90),calibrated=int(0.93<=k/n<=0.97)))
    return rows
def load_reps(path):
    R=json.load(open(path)); reps={}
    for x in R: reps.setdefault(x['rep'],{}).setdefault(x['d'],([0]*5,[0]*5)); reps[x['rep']][x['d']][0][x['ip']]=x['k']; reps[x['rep']][x['d']][1][x['ip']]=x['n']
    return [{d:(np.array(v[0],float),np.array(v[1],float)) for d,v in reps[r].items()} for r in sorted(reps)]
