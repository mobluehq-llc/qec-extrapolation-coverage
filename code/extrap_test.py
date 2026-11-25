# Held-out extrapolation test on the stim data: fit e=A*(p/pth)^((d+1)/2), predict held-out cells, check against measured (binomial) with delta-method CI.
import json, numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
rows=json.load(open('repro_full.json'))
def P_of(th,d,p,r):
    lA,lq=th  # e = exp(lA) * (p*exp(lq))^((d+1)/2)  (q=1/pth)
    e=np.exp(lA+(d+1)/2*(np.log(p)+lq)); e=np.minimum(e,0.499)
    return (1-(1-2*e)**r)/2
def nll(th,R):
    s=0
    for x in R:
        P=P_of(th,x['d'],x['p'],x['r']); n,k=x['shots'],x['errors']
        s-= k*np.log(P)+(n-k)*np.log1p(-P)
    return s
def fit(R):
    best=minimize(nll,[np.log(0.1),np.log(100)],args=(R,),method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-9,maxiter=4000))
    th=best.x; h=1e-4; H=np.zeros((2,2))
    f=lambda t:nll(t,R)
    for i in range(2):
        for j in range(2):
            ei=np.eye(2)[i]*h; ej=np.eye(2)[j]*h
            H[i,j]=(f(th+ei+ej)-f(th+ei-ej)-f(th-ei+ej)+f(th-ei-ej))/(4*h*h)
    return th,np.linalg.inv(H)
def report(name,train,test):
    th,C=fit(train)
    print(f'--- {name}: fit A={np.exp(th[0]):.4f} pth={1/np.exp(th[1]):.4f}  (n_train={len(train)})')
    zs=[]
    for x in test:
        P=P_of(th,x['d'],x['p'],x['r']); g=np.array([(P_of(th+np.eye(2)[i]*1e-5,x['d'],x['p'],x['r'])-P_of(th-np.eye(2)[i]*1e-5,x['d'],x['p'],x['r']))/2e-5 for i in range(2)])
        sd_pred=np.sqrt(g@C@g); n,k=x['shots'],x['errors']; Pm=k/n; sd_m=np.sqrt(max(k,1))/n
        z=(Pm-P)/np.sqrt(sd_pred**2+sd_m**2); zs.append(z)
        print(f"  held-out d={x['d']} p={x['p']:.0e}: pred P_shot={P:.3e} (fit-sd {sd_pred/P*100:.0f}% rel)  measured {Pm:.3e} ({k} errs)  z={z:+.2f}")
    return zs
tr=lambda f:[x for x in rows if f(x)]
report('A: train d=3,5 all p -> predict d=7',tr(lambda x:x['d']<7),tr(lambda x:x['d']==7))
report('B: train p>=3e-3 all d -> predict p=1e-3,2e-3 (low-p extrapolation)',tr(lambda x:x['p']>=3e-3),tr(lambda x:x['p']<3e-3))
report('C: train p>=3e-3, d=3,5 -> predict d=7, p=1e-3,2e-3 (both)',tr(lambda x:x['p']>=3e-3 and x['d']<7),tr(lambda x:x['p']<3e-3 and x['d']==7))
