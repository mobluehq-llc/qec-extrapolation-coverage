# BW-2 ratio + table helpers for the below-threshold amendment (reads frozen outputs; no new fitting choices)
import numpy as np,csv,sys
from coverage_bw import *
from analyze_bw import load_truth
R='../results'
T=load_truth(); truth={k:x['truth'] for k,x in T.items()}
reps=load_reps(f'{R}/fit_bw.json')
rat={}
for data in reps:
    for tp in TARGETS:
        o=fit_all({d:data[d] for d in (3,5,7)},tp)
        for (rc,d),(mu,var) in o.items():
            rat.setdefault((rc,d,tp),[]).append(np.exp(mu)/truth[(d,tp)])
print('median prediction/truth (below-threshold window)')
for rc in ('R1','R2','R3','R4','R5'):
    for d in (3,5,7):
        if (rc,d,1e-3) in rat: print(rc,d,' '.join(f'{np.median(rat[(rc,d,t)]):.3f}' for t in TARGETS))
old={(r['recipe'],int(r['d']),float(r['target'])):float(r['cov']) for r in csv.DictReader(open(f'{R}/L_coverage.csv'))}
new={(r['recipe'],int(r['d']),float(r['target'])):float(r['cov']) for r in csv.DictReader(open(f'{R}/L_bw_coverage.csv'))}
print('R1 strictly higher in',sum(new[k]>old[k] for k in new if k[0]=='R1'),'of 12')
for rc in ('R1','R2','R3','R4','R5'):
    print(rc,'new vs old cov:',' '.join(f'{new[k]:.2f}/{old[k]:.2f}' for k in sorted(new) if k[0]==rc))
for rc,n in (('R2',None),('R3',None)):
    ks=[k for k in new if k[0]==rc]; far=[k for k in ks if k[2]<=1e-4]
    print(rc,'miscal all',sum(new[k]<.9 for k in ks),len(ks),'far',sum(new[k]<.9 for k in far),len(far))
b=[k for k in new if k[0] in('R1','R5')]; print('R1/R5 cells in 93-97:',sum(.93<=new[k]<=.97 for k in b),'of',len(b))
