# one-shot holdout: d=9 at p=1e-4. truth = inverse-variance combination of independent runs of the validated method.
import json,glob,csv,numpy as np
from coverage import *
R='../results'
runs=[json.load(open(f)) for f in sorted(glob.glob(f'{R}/holdout_truth_d9_p1e-4_s*.json'))]
w=np.array([1/r['se']**2 for r in runs]); t=np.array([r['truth'] for r in runs]); truth=float((w*t).sum()/w.sum()); se=float(w.sum()**-0.5)
print('d=9 p=1e-4 truth runs',[(r['truth'],r['se']) for r in runs],'combined',truth,'+-',se,'rel',se/truth)
reps3=load_reps(f'{R}/fit_d357.json'); reps9=load_reps(f'{R}/fit_d9.json')
reps=[{**a,**b} for a,b in zip(reps3,reps9)]
rows=coverage(reps,{(9,1e-4):truth},ds=(3,5,9),targets=[1e-4])
w_=csv.DictWriter(open(f'{R}/L_holdout_d9.csv','w'),fieldnames=list(rows[0])); w_.writeheader(); w_.writerows(rows)
for r in rows: print(r)
get=lambda rc:[r for r in rows if r['recipe']==rc and r['d']==9][0]
r1,r5=get('R1'),get('R5')
print(f"HOLDOUT d=9 p=1e-4: R1 coverage {r1['cov']:.3f} CP[{r1['cp_lo']:.3f},{r1['cp_hi']:.3f}]; R5 coverage {r5['cov']:.3f} CP[{r5['cp_lo']:.3f},{r5['cp_hi']:.3f}] -> "+("CONFIRMED (both <90%)" if r1['cov']<0.9 and r5['cov']<0.9 else "NOT CONFIRMED"))
# direction + magnitude for information
for rc in ('R1','R2','R3','R4','R5'):
    lr=[]
    for data in reps:
        o=fit_all({d:data[d] for d in (3,5,9)},1e-4); lr.append(o[(rc,9)][0]-np.log(truth))
    print(rc,'median pred/truth %.3g'%np.exp(np.median(lr)))
