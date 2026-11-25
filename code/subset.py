# Ground-truth for a SURROGATE noise model: each matching-graph edge of the stim rotated-surface-code memory DEM fires iid with prob p.
# Decoder = pymatching with weights fixed from the circuit-level DEM at p=1e-3 (same decoder at every p).
# p_L(p) = sum_w Binom(w;N,p) f_w, f_w = failure fraction given exactly w random edges fire (stratified subset sampling; exact for iid-uniform edges).
import sys, json, time, numpy as np, stim, pymatching
from scipy.stats import binom
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
def build(d):
    c=stim.Circuit.generated("surface_code:rotated_memory_z",distance=d,rounds=d,after_clifford_depolarization=1e-3,
        before_round_data_depolarization=1e-3,before_measure_flip_probability=1e-3,after_reset_flip_probability=1e-3)
    m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
    D=m.num_detectors; E=[]
    for u,v,a in m.edges():
        E.append((u,-1 if v is None else v, 1 if a.get('fault_ids') else 0))
    return m,D,E
def tables(D,E):
    N=len(E); H=np.zeros((N,D),dtype=np.uint8); L=np.zeros(N,dtype=np.uint8)
    for i,(u,v,o) in enumerate(E):
        H[i,u]^=1
        if v>=0: H[i,v]^=1
        L[i]=o
    return H,L
def fail_frac_subset(m,H,L,w,M):
    N=H.shape[0]; idx=rng.integers(0,N,size=(M,w))
    for _ in range(50):
        s=np.sort(idx,axis=1); bad=(s[:,1:]==s[:,:-1]).any(axis=1) if w>1 else np.zeros(M,bool)
        if not bad.any(): break
        idx[bad]=rng.integers(0,N,size=(bad.sum(),w))
    dets=np.zeros((M,H.shape[1]),dtype=np.uint8); obs=np.zeros(M,dtype=np.uint8)
    for j in range(w): dets^=H[idx[:,j]]; obs^=L[idx[:,j]]
    pred=m.decode_batch(dets); return int((pred[:,0]!=obs).sum()), M
def fail_direct(m,H,L,p,M):
    N=H.shape[0]; fire=rng.random((M,N))<p   # M*N bool; chunked by caller
    dets=(fire.astype(np.uint8)@H)&1; obs=(fire.astype(np.uint8)@L)&1
    pred=m.decode_batch(dets.astype(np.uint8)); return int((pred[:,0]!=obs).sum()), M
if __name__=='__main__':
    out={}
    for d in (3,5,7,9):
        m,D,E=build(d); H,L=tables(D,E); N=len(E); t0=time.time()
        res={'N':N,'D':D,'subset':{},'direct':{}}
        # direct MC at moderate p (chunked)
        for p in (0.004,0.006,0.01,0.015):
            k=n=0; target=300
            while k<target and n<4_000_000:
                kk,nn=fail_direct(m,H,L,p,20000); k+=kk; n+=nn
            res['direct'][str(p)]=[k,n]
        # subset sampling w=1..12
        for w in range(1,13):
            k=n=0
            while n<400_000 and (k<400 or w<=(d+1)//2+1):
                kk,nn=fail_frac_subset(m,H,L,w,50000); k+=kk; n+=nn
                if k>=3000: break
            res['subset'][str(w)]=[k,n]
        res['secs']=time.time()-t0; out[str(d)]=res
        print(d,N,D,'secs',round(res['secs']),flush=True)
        json.dump(out,open('subset_results.json','w'))
