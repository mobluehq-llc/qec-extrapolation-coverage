# Control (b): edge-iid surrogate from F2 (positive control with known behaviour). Replicate window data by direct simulation of the
# surrogate; truth from F2's stratified subset results (subset_results.json + boost7.json) at targets where truncation tail<1e-3.
import sys,json,numpy as np
from multiprocessing import Pool
sys.argv=[sys.argv[0]]+sys.argv[1:]
import subset as S
PS=[4e-3,6e-3,8e-3,1e-2,1.5e-2]; _c={}
def get(d):
    if d not in _c:
        m,D,E=S.build(d); H,L=S.tables(D,E); _c[d]=(m,H,L)
    return _c[d]
def one(a):
    rep,d,ip=a; m,H,L=get(d); ss=np.random.SeedSequence([12,rep,d,ip]); S.rng=np.random.default_rng(ss); p=PS[ip]; k=n=0
    while k<300:
        kk,nn=S.fail_direct(m,H,L,p,5000 if d<7 else 2000); k+=kk; n+=nn
    return dict(rep=rep,d=d,ip=ip,k=k,n=n)
if __name__=='__main__':
    nproc=int(sys.argv[1]); nrep=int(sys.argv[2]); out=sys.argv[3]
    tasks=[(r,d,ip) for r in range(nrep) for d in (3,5,7) for ip in range(5)]
    with Pool(nproc) as P: res=P.map(one,tasks,chunksize=5)
    json.dump(res,open(out,'w')); print('done',len(res))
