# Control (a): DATA-ONLY noise (before_round_data_depolarization=p, perfect gates/measurement). One window dataset per (d,ip), >=300 errors or maxshots cap.
import sys,json,time,numpy as np,stim,pymatching
from multiprocessing import Pool
from truth import circuit
PS=[4e-3,6e-3,8e-3,1e-2,1.5e-2]
def one(a):
    d,ip,cap=a; c=circuit(d,PS[ip],'dataonly'); m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
    s=c.compile_detector_sampler(seed=int(np.random.SeedSequence([13,d,ip]).generate_state(1,dtype=np.uint64)[0]%2**63)); k=n=0; B=200000
    while k<300 and n<cap:
        x,o=s.sample(B,separate_observables=True); k+=int((m.decode_batch(x)[:,0]!=o[:,0]).sum()); n+=B
    return dict(d=d,ip=ip,k=k,n=n)
if __name__=='__main__':
    nproc=int(sys.argv[1]); out=sys.argv[2]; cap=float(sys.argv[3])
    tasks=[(d,ip,cap) for d in (7,5,3) for ip in range(5)]
    with Pool(nproc) as P: res=P.map(one,tasks,chunksize=1)
    json.dump(res,open(out,'w')); print(res)
