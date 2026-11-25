"""Item 8(b): R2E/R3E = R2/R3 fitted to P-hat/E(p), E = exact zero-fault probability from the undecomposed DEM.
Run: python code/zerofault_r2e.py  (from the L folder or code/). Needs stim, numpy, scipy. Writes results/L_zerofault.json, L_zerofault_coverage.csv."""
import os,sys,json,csv,numpy as np,stim
os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,'.')
def circuit(d,p,kind='full'):  # verbatim copy of truth.circuit (truth.py imports pymatching, absent on a compute host)
    return stim.Circuit.generated("surface_code:rotated_memory_z",distance=d,rounds=d,
        after_clifford_depolarization=p,before_measure_flip_probability=p,after_reset_flip_probability=p)
import coverage as A, coverage_bw as B, analyze_bw as AB
R='../results'; N=300; B90=270; TG=A.TARGETS
def logE(d,p):
    dem=circuit(d,p,'full').detector_error_model(); q=np.array([i.args_copy()[0] for i in dem.flattened() if i.type=='error'])
    return float(np.log1p(-q).sum())
def fit(data,tp,LE,LEt,PS,scaled=False):
    """data[d]=(k[5],n[5]); LE[d]=ln E at window pts; LEt[d]=ln E at target; returns {(rec,d):(mu,var)} for ln P at target"""
    X0=np.log(PS); xs=np.log(tp); out={}; Y={};W={}
    for d in data:
        k,n=data[d]; P=k/n; Y[d]=np.log(P)-LE[d]; W[d]=k/(1-P)
    def sc(C,c2,dof): return C*max(1.0,c2/dof) if scaled and dof>0 else C
    for d in data:
        h=(d+1)/2; b,C,c2,dof=A.wls(np.ones((5,1)),Y[d]-h*X0,W[d]); C=sc(C,c2,dof); out[('R2E',d)]=(b[0]+h*xs+LEt[d],C[0,0])
    Xs=[];ys=[];ws=[]
    for d in data:
        h=(d+1)/2; Xs.append(np.c_[np.ones(5),np.full(5,h)]); ys.append(Y[d]-h*X0); ws.append(W[d])
    b,C,c2,dof=A.wls(np.vstack(Xs),np.concatenate(ys),np.concatenate(ws)); C=sc(C,c2,dof)
    for d in data:
        h=(d+1)/2; g=np.array([1,h]); out[('R3E',d)]=(g@b+h*xs+LEt[d],g@C@g)
    return out
T=AB.load_truth(); truth={k:v['truth'] for k,v in T.items()}
frozen={'w1':{(r['recipe'],int(r['d']),float(r['target'])):int(r['hits']) for r in csv.DictReader(open(f'{R}/L_coverage.csv'))},
        'w2':{(r['recipe'],int(r['d']),float(r['target'])):int(r['hits']) for r in csv.DictReader(open(f'{R}/L_bw_coverage.csv'))}}
LET={(d,tp):logE(d,tp) for d in (3,5,7) for tp in TG}
out={'logE_targets':{f'{d}|{tp:g}':v for (d,tp),v in LET.items()}}; rows=[]
for w,M,f in (('w1',A,'fit_d357.json'),('w2',B,'fit_bw.json')):
    PS=M.PS; LE={d:np.array([logE(d,p) for p in PS]) for d in (3,5,7)}
    out[f'logE_window_{w}']={str(d):list(v) for d,v in LE.items()}
    reps=M.load_reps(f'{R}/{f}'); assert len(reps)==300
    for scaled in (False,True):
        # sanity: plain R2/R3 reproduce frozen hits (primary only; frozen scaled csv checked too)
        if not scaled:
            for tp in TG:
                for rc in ('R2','R3'):
                    for d in (3,5,7):
                        h=sum(abs(np.log(truth[(d,tp)])-M.fit_all({dd:x[dd] for dd in (3,5,7)},tp)[(rc,d)][0])<=1.96*np.sqrt(M.fit_all({dd:x[dd] for dd in (3,5,7)},tp)[(rc,d)][1]) for x in reps)
                        assert h==frozen[w][(rc,d,tp)],(w,rc,d,tp,h,frozen[w][(rc,d,tp)])
        for tp in TG:
            LEt={d:LET[(d,tp)] for d in (3,5,7)}
            res=[fit({dd:x[dd] for dd in (3,5,7)},tp,LE,LEt,PS,scaled) for x in reps]
            for rc in ('R2E','R3E'):
                for d in (3,5,7):
                    mu=np.array([r[(rc,d)][0] for r in res]); sd=1.96*np.sqrt([r[(rc,d)][1] for r in res]); t0=truth[(d,tp)]; se=T[(d,tp)]['se']
                    hit=lambda t:(np.abs(np.log(t)-mu)<=sd).sum()
                    k=int(hit(t0)); gk=[int(hit(t)) for t in np.linspace(t0-1.96*se,t0+1.96*se,81)]
                    lab='below' if max(gk)<B90 else ('above' if min(gk)>=B90 else 'undecided')
                    rows.append(dict(window=w,scaled=int(scaled),recipe=rc,d=d,target=tp,hits=k,cov=k/N,kmin=min(gk),kmax=max(gk),label=lab,
                        median_ratio=float(np.median(np.exp(mu)/t0)),median_halfwidth_ln=float(np.median(sd)),point_below90=int(k<B90)))
w=csv.DictWriter(open(f'{R}/L_zerofault_coverage.csv','w'),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def sel(w,rc,scaled=0,ds=(3,5)): return [r for r in rows if r['window']==w and r['recipe']==rc and r['scaled']==scaled and r['d'] in ds]
P1=[r for rc in ('R2E','R3E') for r in sel('w1',rc)]; p1ok=all(0.90<=r['median_ratio']<=1.10 for r in P1)
P2n=sum(r['point_below90'] for r in P1); P3=sel('w2','R2E'); P3n=sum(r['label']=='below' for r in P3)
verd=dict(P1=dict(cells=len(P1),in_band=sum(0.90<=r['median_ratio']<=1.10 for r in P1),ratio_range=[min(r['median_ratio'] for r in P1),max(r['median_ratio'] for r in P1)],MET=bool(p1ok)),
          P2=dict(below90=P2n,cells=len(P1),MET=bool(P2n<=8)),P3=dict(robust_below=P3n,cells=len(P3),MET=bool(P3n<=2)))
out['verdicts']=verd; json.dump(out,open(f'{R}/L_zerofault.json','w'),indent=1)
print(json.dumps(verd,indent=1))
for r in rows:
    if r['scaled']==0: print(r['window'],r['recipe'],r['d'],f"{r['target']:g}",r['hits'],f"{r['median_ratio']:.3f}",r['kmin'],r['kmax'],r['label'],f"hw={r['median_halfwidth_ln']:.3f}")
