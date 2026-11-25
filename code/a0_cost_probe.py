import sys, time, math, numpy as np, stim, hashlib
sys.path.insert(0, "code")
from truth import Model, circuit
c=circuit(9,1e-4)
loc=0
for ins in c.flattened():
    n=ins.name; k=len(ins.target_groups()) if hasattr(ins,'target_groups') else len(ins.targets_copy())
    if n in("X_ERROR","Z_ERROR"): loc+=k
    elif n=="DEPOLARIZE1": loc+=3*k
    elif n=="DEPOLARIZE2": loc+=15*k
print("pauli_components",loc)
M=Model(9,1e-4,1); N=M.N
print("N",N,"unique det sets",len(set(M.dets)))
rng=np.random.default_rng(7)
for w in (2,3,4):
    B=20000
    idx=np.sort(rng.integers(0,N,size=(B,w)),axis=1); t=time.time(); f=M.fail(idx); dt=time.time()-t
    print("random w",w,"us/decode",1e6*dt/B,"fails",int(f.sum()))
    idx,_=(M.sample_walk(B,w)[0],None) if False else (None,None)
    try:
        r=M.sample_walk(B,w); idx=r[0] if isinstance(r,tuple) else r
        t=time.time(); f=M.fail(idx); dt=time.time()-t
        print("walk(local) w",w,"us/decode",1e6*dt/B,"fails",int(f.sum()))
    except Exception as e: print("walk err",e)
print("sha_circuit",hashlib.sha256(str(c).encode()).hexdigest()[:16])
