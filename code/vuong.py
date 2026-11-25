# EXPLORATORY: Vuong (1989) non-nested comparison of R1..R4 pooled over d=3,5,7, Gaussian-in-log-P likelihood with known sigma, BIC-corrected.
import numpy as np,itertools,json,sys
from coverage import PS,X0,load_reps,wls
def lls(data):
    out={}
    res={m:[] for m in ('R1','R2','R3','R4')}; npar={'R1':6,'R2':3,'R3':2,'R4':9}
    Y={};W={}
    for d,(k,n) in data.items(): P=k/n; Y[d]=np.log(P); W[d]=k/(1-P)
    def ll(y,yh,w): return -0.5*(w*(y-yh)**2-np.log(w)+np.log(2*np.pi))
    for d in data:
        y,w=Y[d],W[d]
        X=np.c_[np.ones(5),X0]; b=wls(X,y,w)[0]; res['R1'].append(ll(y,X@b,w))
        h=(d+1)/2; b=wls(np.ones((5,1)),y-h*X0,w)[0]; res['R2'].append(ll(y,b[0]+h*X0,w))
        X=np.c_[np.ones(5),PS/1e-2,(PS/1e-2)**2]; b=wls(X,y-(d/2)*X0,w)[0]; res['R4'].append(ll(y,X@b+(d/2)*X0,w))
    Xs=[];ys=[];ws=[]
    for d in data: h=(d+1)/2; Xs.append(np.c_[np.ones(5),np.full(5,h)]); ys.append(Y[d]-h*X0); ws.append(W[d])
    X=np.vstack(Xs);y=np.concatenate(ys);w=np.concatenate(ws); b=wls(X,y,w)[0]; fit=X@b
    res['R3'].append(ll(y,fit,w))
    return {m:np.concatenate(v) for m,v in res.items()},npar
def run(path,out):
    reps=load_reps(path); tally={}
    for data in reps:
        L,npar=lls({d:data[d] for d in (3,5,7)}); n=len(L['R1'])
        for a,b in itertools.combinations(L,2):
            m=L[a]-L[b]-(npar[a]-npar[b])/2*np.log(n)/n
            z=np.sqrt(n)*m.mean()/(m.std(ddof=1)+1e-12)
            v='undecided' if abs(z)<1.96 else (a if z>0 else b)
            tally.setdefault(f'{a}_vs_{b}',{}).setdefault(v,0); tally[f'{a}_vs_{b}'][v]+=1
    json.dump(tally,open(out,'w'),indent=1)
    for k,v in tally.items(): print(k,v)
if __name__=='__main__': run(sys.argv[1],sys.argv[2])
