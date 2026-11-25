import json,sys,numpy as np,csv
from scipy.stats import binom
from coverage import *
R='../results'
def dataset(rows):
    d={}
    for x in rows: d.setdefault(x['d'],([0]*5,[0]*5)); d[x['d']][0][x['ip']]=x['k']; d[x['d']][1][x['ip']]=x['n']
    return {k:(np.array(v[0],float),np.array(v[1],float)) for k,v in d.items()}
def control_a():
    reps=load_reps(f'{R}/fit_d357.json'); DO=dataset(json.load(open(f'{R}/ctrl_a_dataonly.json')))
    res={}  # (recipe,d,ip,kind)->[rej,total]
    for r in range(150):
        fit=reps[r]; fresh=reps[r+150]
        for ip,p in enumerate(PS):
            o=fit_all(fit,p)
            for (rc,d),(mu,var) in o.items():
                if rc=='R5': continue
                for kind,D in (('dataonly',DO),('circuit_fresh',fresh)):
                    k,n=D[d][0][ip],D[d][1][ip]; P=k/n; z=(np.log(P)-mu)/np.sqrt(var+(1-P)/k)
                    c=res.setdefault((rc,d,ip,kind),[0,0]); c[0]+=int(abs(z)>1.96); c[1]+=1
    rows=[dict(recipe=a,d=b,window_p=PS[c],kind=k,rejected=v[0],tests=v[1],reject_rate=v[0]/v[1]) for (a,b,c,k),v in sorted(res.items())]
    w=csv.DictWriter(open(f'{R}/L_control_dataonly.csv','w'),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    out=[]
    for rc in ('R1','R2','R3','R4'):
        for kind in ('dataonly','circuit_fresh'):
            s=[r for r in rows if r['recipe']==rc and r['kind']==kind]
            out.append(f"CONTROL-A {rc} {kind}: reject rate {sum(r['rejected'] for r in s)/sum(r['tests'] for r in s):.3f}; cells with reject rate>=0.95: {sum(r['reject_rate']>=0.95 for r in s)}/{len(s)}")
    out.append('CONTROL-A data-only measured P/ circuit-fit: see csv; data-only points (d,p,k,n): '+str({d:(DO[d][0].tolist(),DO[d][1].tolist()) for d in DO}))
    return out
def surrogate_truth():
    S=json.load(open(f'{R}/F2_subset_results.json'));B=json.load(open(f'{R}/F2_boost7.json'))
    for w,(k,n) in B.items():
        k0,n0=S['7']['subset'][w]; S['7']['subset'][w]=[k0+k,n0+n]
    T={};tails={}
    for d in S:
        r=S[d]; N=r['N']
        for tp in TARGETS:
            wmin=min(int(w) for w,(k,n) in r['subset'].items() if k>0); mu=0;var=0
            for w,(k,n) in r['subset'].items():
                w=int(w)
                if w<wmin: continue
                pm=binom.pmf(w,N,tp); f=k/n; mu+=pm*f; var+=pm**2*max(f*(1-f),0)/n
            tail=binom.sf(max(int(w) for w in r['subset']),N,tp)
            tails[(int(d),tp)]=tail
            if tail<1e-3 and mu>0: T[(int(d),tp)]=mu
    return T,tails
def control_b():
    T,tails=surrogate_truth(); reps=load_reps(f'{R}/ctrl_b_surrogate.json')
    rows=coverage(reps,T); w=csv.DictWriter(open(f'{R}/L_control_surrogate.csv','w'),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    out=[f"CONTROL-B surrogate truth cells used (tail<1e-3): {sorted(T)}; replicates {len(reps)}"]
    for rc in ('R1','R2','R3','R4','R5'):
        s=[r for r in rows if r['recipe']==rc]
        out.append(f"CONTROL-B {rc}: miscalibrated {sum(r['miscal'] for r in s)}/{len(s)} cells; median coverage {np.median([r['cov'] for r in s]):.3f}")
    far=[r for r in rows if r['recipe']=='R1' and r['d'] in (5,7) and r['target']<=1e-5]
    out.append(f"CONTROL-B R1 d=5,7 targets<=1e-5: miscalibrated {sum(r['miscal'] for r in far)}/{len(far)}")
    return out
if __name__=='__main__':
    which=sys.argv[1]; out=control_a() if which=='a' else control_b(); open(f'{R}/L_control_{which}.txt','w').write('\n'.join(out)); print('\n'.join(out))
