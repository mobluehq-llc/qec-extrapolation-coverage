# Amendment 2: rolling-origin discrepancy inflation for R2/R3 (and R1 side check). Stored data only.
import json,csv,numpy as np
import coverage_bw as CB, coverage as CS
from analyze_bw import load_truth
R='../results'; TARGETS=[1e-3,1e-4,1e-5,1e-6]; DS=(3,5,7); Z=1.96
def fit(data,PS,idx,tp,scaled=False):
    """fit R1,R2,R3 on points idx (list); return {(rc,d):(mu,var)} at tp and 'pred(p)' function via closure-free recompute"""
    X0=np.log(PS[idx]); xs=np.log(tp); out={}; Y={};W={}
    for d in DS:
        k,n=data[d]; P=k[idx]/n[idx]; Y[d]=np.log(P); W[d]=k[idx]/(1-P)
    sc=lambda C,c2,dof: C*max(1.0,c2/dof) if scaled and dof>0 else C
    m=len(idx)
    for d in DS:
        X=np.c_[np.ones(m),X0]; b,C,c2,dof=CB.wls(X,Y[d],W[d]); C=sc(C,c2,dof); g=np.array([1,xs]); out[('R1',d)]=(g@b,g@C@g)
        h=(d+1)/2; b,C,c2,dof=CB.wls(np.ones((m,1)),Y[d]-h*X0,W[d]); C=sc(C,c2,dof); out[('R2',d)]=(b[0]+h*xs,C[0,0])
    Xs=[];ys=[];ws=[]
    for d in DS:
        h=(d+1)/2; Xs.append(np.c_[np.ones(m),np.full(m,h)]); ys.append(Y[d]-h*X0); ws.append(W[d])
    b,C,c2,dof=CB.wls(np.vstack(Xs),np.concatenate(ys),np.concatenate(ws)); C=sc(C,c2,dof)
    for d in DS:
        h=(d+1)/2; g=np.array([1,h]); out[('R3',d)]=(g@b+h*xs,g@C@g)
    return out
def run(reps,PS,truth,leave,scaled=False):
    hold=len(PS)-1 if leave=='high' else 0
    idx=[i for i in range(len(PS)) if i!=hold]; allidx=list(range(len(PS)))
    cnt={};cnt0={};wl=[];wlin={}
    for data in reps:
        e={}
        sub=fit(data,PS,idx,PS[hold],scaled)
        for (rc,d),(mu,var) in sub.items():
            k,n=data[d]; e[(rc,d)]=np.log(k[hold]/n[hold])-mu
        for tp in TARGETS:
            full=fit(data,PS,allidx,tp,scaled)
            for (rc,d),(mu,var) in full.items():
                t=truth.get((d,tp))
                if t is None: continue
                s0=np.sqrt(var); s1=np.sqrt(var+e[(rc,d)]**2); lt=np.log(t)
                key=(rc,d,tp)
                cnt[key]=cnt.get(key,0)+int(abs(lt-mu)<=Z*s1); cnt0[key]=cnt0.get(key,0)+int(abs(lt-mu)<=Z*s0)
                wl.append((rc,s1/s0)); 
                lw=lambda s: np.exp(mu+Z*s)-np.exp(mu-Z*s)
                wlin.setdefault(rc,[]).append(lw(s1)/lw(s0))
    rows=[]
    for key in sorted(cnt):
        k=cnt[key];n=len(reps);lo,hi=CB.cp(k,n);k0=cnt0[key]
        rows.append(dict(recipe=key[0],d=key[1],target=key[2],hits=k,reps=n,cov=k/n,cp_lo=lo,cp_hi=hi,cov_uninflated=k0/n,hits_uninflated=k0))
    wd={rc:(float(np.median([w for r,w in wl if r==rc])),float(np.median(wlin[rc]))) for rc in wlin}
    return rows,wd
def load(path,PS):
    return CB.load_reps(path)
def summarize(rows,wd,label):
    out=[]
    for rc in ('R2','R3','R1'):
        c=[r for r in rows if r['recipe']==rc]
        band=sum(.93<=r['cov']<=.97 for r in c); lo=sum(r['cov']<.90 for r in c); hi=sum(r['cov']>.97 for r in c)
        b0=sum(.93<=r['cov_uninflated']<=.97 for r in c); l0=sum(r['cov_uninflated']<.90 for r in c)
        out.append(f"{label} {rc}: inflated in-band(93-97) {band}/{len(c)}, <90 {lo}, >97 {hi}, median cov {np.median([r['cov'] for r in c]):.3f}; uninflated in-band {b0}, <90 {l0}, median cov {np.median([r['cov_uninflated'] for r in c]):.3f}; median width ratio lnscale {wd[rc][0]:.2f}x linear {wd[rc][1]:.2f}x")
    c=[r for r in rows if r['recipe'] in ('R2','R3')]
    out.append(f"{label} BAR cells R2+R3 in band: {sum(.93<=r['cov']<=.97 for r in c)}/{len(c)} (need >=18)")
    far=[r for r in c if r['target']<=1e-4]
    out.append(f"{label} far targets (<=1e-4) R2+R3 in band: {sum(.93<=r['cov']<=.97 for r in far)}/{len(far)}, <90: {sum(r['cov']<.90 for r in far)}")
    return out
if __name__=='__main__':
    T=load_truth(); truth={(d,p):x['truth'] for (d,p),x in T.items()}
    W={'BW':(CB.PS,f'{R}/fit_bw.json'),'STRADDLE':(CS.PS,f'{R}/fit_d357.json')}
    lines=[];allrows=[];bar={}
    for name,(PS,path) in W.items():
        reps=CB.load_reps(path); assert len(reps)==300
        for leave,scaled,tag in (('high',False,'CONFIRMATORY'),('high',True,'EXPLORATORY-scaled'),('low',False,'EXPLORATORY-leave-lowest')):
            rows,wd=run(reps,PS,truth,leave,scaled)
            for r in rows: r.update(window=name,variant=tag)
            allrows+=rows; lines+=summarize(rows,wd,f"[{name} {tag}]")
            if tag=='CONFIRMATORY':
                c=[r for r in rows if r['recipe'] in('R2','R3')]; bar[name]=sum(.93<=r['cov']<=.97 for r in c)
                # sanity: uninflated primary must equal stored coverage tables
                st=list(csv.DictReader(open(f'{R}/L_bw_coverage.csv' if name=='BW' else f'{R}/L_coverage.csv')))
                sd={(x['recipe'],int(x['d']),float(x['target'])):int(x['hits']) for x in st}
                bad=[k for k in ((r['recipe'],r['d'],r['target']) for r in rows if r['recipe'] in('R1','R2','R3')) if False]
                mism=[(r['recipe'],r['d'],r['target'],r['hits_uninflated'],sd.get((r['recipe'],r['d'],r['target']))) for r in rows if sd.get((r['recipe'],r['d'],r['target']))!=r['hits_uninflated']]
                lines.append(f"[{name}] sanity uninflated vs stored table mismatches: {len(mism)} {mism[:3]}")
    lines.append(f"VERDICT (frozen rule, both windows, need >=18/24 each): BW {bar['BW']}/24, STRADDLE {bar['STRADDLE']}/24 -> {'PASS' if min(bar.values())>=18 else 'FAIL'}")
    w=csv.DictWriter(open(f'{R}/L_inflate_coverage.csv','w'),fieldnames=list(allrows[0])); w.writeheader(); w.writerows(allrows)
    open(f'{R}/L_inflate_verdicts.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
