import json,glob,csv,sys,numpy as np
from coverage_bw import *
R='../results'
def load_truth():
    T={}
    for f in glob.glob(f'{R}/truth_d*_p*.json'):
        x=json.load(open(f)); T[(x['d'],x['p'])]=x
    for f in glob.glob(f'{R}/truthextra_d*_p*.json'):   # independent second run of the same cell: inverse-variance combine
        y=json.load(open(f)); x=T[(y['d'],y['p'])]; w1,w2=1/x['se']**2,1/y['se']**2
        x['truth_run1']=x['truth']; x['se_run1']=x['se']; x['truth_run2']=y['truth']; x['se_run2']=y['se']
        x['truth']=(w1*x['truth']+w2*y['truth'])/(w1+w2); x['se']=(w1+w2)**-0.5
    return T
def write(rows,path):
    w=csv.DictWriter(open(path,'w'),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def verdicts(rows,label):
    out=[]
    def cells(rc,decs): return [r for r in rows if r['recipe']==rc and r['target'] in decs]
    far=[1e-4,1e-5,1e-6]
    for rc in ('R1','R5'):
        c=cells(rc,far); mis=sum(r['miscal'] for r in c); firm=sum(r['miscal_firm'] for r in c)
        out.append(f"{label} H1 component {rc}: {mis}/{len(c)} cells miscalibrated (cov<90%) at >=2 decades (CP-firm: {firm}/{len(c)}); need >= {2*len(c)/3:.2f} -> {'MET' if mis>=2*len(c)/3 else 'NOT MET'}")
    r1=cells('R1',far);r5=cells('R5',far)
    ok=(sum(r['miscal'] for r in r1)>=2*len(r1)/3) and (sum(r['miscal'] for r in r5)>=2*len(r5)/3)
    out.append(f"{label} H1 (R1 and R5 miscalibrated in >=2/3 of cells at >=2 decades): {'SUPPORTED' if ok else 'NOT SUPPORTED'}")
    r1all=[r for r in rows if r['recipe']=='R1']; r2=[r for r in rows if r['recipe']=='R5']
    out.append(f"{label} Claim1 R1 coverage<90% at 2-4 decades: {sum(r['miscal'] for r in r1)}/{len(r1)} cells")
    out.append(f"{label} Claim2 R5 coverage<90%: {sum(r['miscal'] for r in r5)}/{len(r5)} cells (all targets: {sum(r['miscal'] for r in cells('R5',TARGETS))}/{len(cells('R5',TARGETS))})")
    for rc in ('R1','R2','R3','R4','R5'):
        c=[r for r in rows if r['recipe']==rc]; cal=[r for r in c if r['calibrated']]
        out.append(f"{label} Claim3 {rc}: {len(cal)}/{len(c)} cells with coverage in [93%,97%]; median coverage {np.median([r['cov'] for r in c]):.3f}; miscalibrated {sum(r['miscal'] for r in c)}/{len(c)}")
    return out
if __name__=='__main__':
    T=load_truth(); truth={(d,p):x['truth'] for (d,p),x in T.items()}
    print('truth cells',len(truth))
    reps=load_reps(f'{R}/fit_bw.json'); print('replicates',len(reps))
    lines=[]
    for scaled,name in ((False,'primary(absolute sigma)'),(True,'sensitivity(chi2/dof-scaled)')):
        rows=coverage(reps,truth,scaled=scaled)
        write(rows,f'{R}/L_bw_coverage{"_scaled" if scaled else ""}.csv')
        lines+=verdicts(rows,name)
    open(f'{R}/L_bw_verdicts.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
