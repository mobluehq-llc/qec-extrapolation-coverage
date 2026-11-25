# replicate window datasets: for each (rep,d,ip) simulate stim circuit + pymatching(weights at that p) until >=300 logical errors.
# seed = SeedSequence([14,rep,d,ip]); chunk size fixed per (d,ip) from a pilot (chunks.json) so that <~20 errors overshoot.
import sys,json,time,numpy as np,stim,pymatching
from multiprocessing import Pool
from truth import circuit
PS=[1.5e-3,2e-3,2.5e-3,3e-3,4e-3]   # BELOW-THRESHOLD window (PLAN.md Below-threshold amendment)
_cache={}
def get(d,ip,kind):
    k=(d,ip,kind)
    if k not in _cache:
        c=circuit(d,PS[ip],kind); m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True)); _cache[k]=(c,m)
    return _cache[k]
def one(a):
    rep,d,ip,chunk,kind=a
    c,m=get(d,ip,kind); ss=np.random.SeedSequence([14,rep,d,ip]); seed=int(ss.generate_state(1,dtype=np.uint64)[0]%(2**63))
    s=c.compile_detector_sampler(seed=seed); k=n=0
    while k<300:
        x,o=s.sample(chunk,separate_observables=True); k+=int((m.decode_batch(x)[:,0]!=o[:,0]).sum()); n+=chunk
    return dict(rep=rep,d=d,ip=ip,k=k,n=n)
if __name__=='__main__':
    mode=sys.argv[1]; out=sys.argv[2]; nproc=int(sys.argv[3]); kind=sys.argv[4] if len(sys.argv)>4 else 'full'
    if mode=='chunks':
        res={}
        for d in (3,5,7):
            for ip,p in enumerate(PS):
                c,m=get(d,ip,kind); s=c.compile_detector_sampler(seed=999); n=1000000 if d<7 else 400000
                x,o=s.sample(n,separate_observables=True); pl=max(1,(m.decode_batch(x)[:,0]!=o[:,0]).sum())/n
                res[f'{d},{ip}']=max(100,int(20/pl))
        json.dump(res,open(out,'w')); print(res)
    else:
        ds=[int(x) for x in mode.split(',')]; nrep=int(sys.argv[5]) if len(sys.argv)>5 else 300
        ch=json.load(open('../results/chunks_bw.json'))
        tasks=[(r,d,ip,ch[f'{d},{ip}'],kind) for r in range(nrep) for d in ds for ip in range(5)]
        t=time.time()
        with Pool(nproc) as P: res=P.map(one,tasks,chunksize=5)
        json.dump(res,open(out,'w')); print('done',len(res),round(time.time()-t))
