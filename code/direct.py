# direct Monte Carlo on the real circuit (stim sampler + pymatching at same p); stop at nerr errors or maxshots
import sys,json,time,stim,pymatching,numpy as np
from truth import circuit
d=int(sys.argv[1]); p=float(sys.argv[2]); seed=int(sys.argv[3]); out=sys.argv[4]; nerr=int(sys.argv[5]); maxshots=float(sys.argv[6]); kind=sys.argv[7] if len(sys.argv)>7 else 'full'
c=circuit(d,p,kind); m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
s=c.compile_detector_sampler(seed=seed); k=n=0; t=time.time(); B=100000 if d<7 else 50000
while k<nerr and n<maxshots:
    x,o=s.sample(B,separate_observables=True); k+=int((m.decode_batch(x)[:,0]!=o[:,0]).sum()); n+=B
json.dump(dict(d=d,p=p,k=k,n=n,secs=time.time()-t,kind=kind),open(out,'w')); print('DIRECT',d,p,k,n,k/n,round(time.time()-t))
