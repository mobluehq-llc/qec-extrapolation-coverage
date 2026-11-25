# L truth: stratified importance sampling on the real circuit DEM (undecomposed mechanisms), pymatching decoder
# with weights from the DEM at the SAME p. p_L(p)=sum_w P(W=w) f_w ; f_w estimated with mixture IS:
#   q = alpha*walk(local growth over detector neighbourhood) + (1-alpha)*Plackett-Luce(prop. to odds r_i)
#   target over ORDERED sequences = prod r_i /(e_w(r) w!)  ;  estimator fail*target/q  (unbiased; W<=~1/(1-alpha))
import sys, json, time, math, numpy as np, stim, pymatching
from scipy.special import gammaln
ALPHA=0.7
def circuit(d,p,kind='full'):
    if kind=='full':
        return stim.Circuit.generated("surface_code:rotated_memory_z",distance=d,rounds=d,
            after_clifford_depolarization=p,before_measure_flip_probability=p,after_reset_flip_probability=p)
    return stim.Circuit.generated("surface_code:rotated_memory_z",distance=d,rounds=d,before_round_data_depolarization=p)
def parse_dem(c):
    dem=c.detector_error_model(); q=[];dets=[];obs=[]
    for ins in dem.flattened():
        if ins.type!='error': continue
        q.append(ins.args_copy()[0]); ds=[t.val for t in ins.targets_copy() if t.is_relative_detector_id()]
        ob=sum(1 for t in ins.targets_copy() if t.is_logical_observable_id())%2
        dets.append(tuple(sorted(ds))); obs.append(ob)
    return np.array(q),dets,np.array(obs,dtype=np.uint8),dem.num_detectors
class Model:
    def __init__(s,d,p,seed,kind='full'):
        s.d,s.p=d,p; s.rng=np.random.default_rng(seed)
        c=circuit(d,p,kind); q,dets,obs,D=parse_dem(c)
        s.q=q; s.N=len(q); s.D=D; s.dets=dets; s.obs=obs
        K=max(len(x) for x in dets); s.K=K
        Hp=np.full((s.N,K),D,dtype=np.int32)
        for i,x in enumerate(dets): Hp[i,:len(x)]=x
        s.Hp=Hp
        s.r=q/(1-q); s.logr=np.log(s.r); s.R=s.r.sum()
        s.m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
        # detector incidence, adjacency (hop<=1, incl self)
        inc=[[] for _ in range(D)]
        for i,x in enumerate(dets):
            for t in x: inc[t].append(i)
        Kd=max(len(x) for x in inc); s.Kd=Kd
        s.deg=np.array([len(x) for x in inc]+[1])
        s.inc=np.zeros((D+1,Kd),dtype=np.int64); s.cumr=np.zeros((D+1,Kd))
        for t in range(D):
            L=inc[t]; s.inc[t,:len(L)]=L; cr=np.cumsum(s.r[L]); s.cumr[t,:len(L)]=cr; s.cumr[t,len(L):]=cr[-1]
        s.cumr[D,:]=1.0
        s.rsum=s.cumr[:,-1].copy()
        Adj=np.zeros((D+1,D+1),dtype=bool)
        for x in dets:
            for a in x:
                for b in x: Adj[a,b]=True
        for a in range(D): Adj[a,a]=True
        s.Adj=Adj
        # stratum probabilities
        s.maxw=min(120,s.N)
        e=np.zeros(s.maxw+1); e[0]=1.0
        for ri in s.r:
            e[1:]=e[1:]+ri*e[:-1]
        s.loge=np.log(np.maximum(e,1e-300)); s.logP0=np.log1p(-q).sum()
        s.logPw=s.loge+s.logP0; s.Pw=np.exp(s.logPw)
        s.lam=q.sum()
    def tail(s,w): return s.Pw[w+1:].sum()
    # ---- samplers
    def sample_pl(s,M,w):
        out=np.empty((M,w),dtype=np.int64); ch=max(1,int(2.5e7//s.N))
        for a in range(0,M,ch):
            m=min(ch,M-a); keys=s.logr[None,:]+s.rng.gumbel(size=(m,s.N))
            idx=np.argpartition(-keys,w-1,axis=1)[:,:w]; kk=np.take_along_axis(keys,idx,1)
            o=np.argsort(-kk,axis=1); out[a:a+m]=np.take_along_axis(idx,o,1)
        return out
    def _addDn(s,Dn,j):
        for t in range(s.K): Dn|=s.Adj[s.Hp[j,t]]
    def sample_walk(s,M,w):
        idx=np.empty((M,w),dtype=np.int64); cr=np.cumsum(s.r)/s.R
        idx[:,0]=np.minimum(np.searchsorted(cr,s.rng.random(M)),s.N-1)
        Dn=np.zeros((M,s.D+1),dtype=bool); s._addDn(Dn,idx[:,0]); ar=np.arange(M)
        for k in range(1,w):
            cs=np.cumsum(Dn[:,:s.D],axis=1); cnt=cs[:,-1]
            u=np.floor(s.rng.random(M)*cnt).astype(np.int64)+1
            dl=(cs>=u[:,None]).argmax(1)
            v=s.rng.random(M)*s.rsum[dl]; pos=(s.cumr[dl]<v[:,None]).sum(1); pos=np.minimum(pos,s.deg[dl]-1)
            j=s.inc[dl,pos]; idx[:,k]=j; s._addDn(Dn,j)
        return idx
    def logq_walk(s,idx):
        M,w=idx.shape; ar=np.arange(M); cr=s.r.sum()
        lq=np.log(s.r[idx[:,0]]/s.R); Dn=np.zeros((M,s.D+1),dtype=bool); s._addDn(Dn,idx[:,0]); bad=np.zeros(M,bool)
        for k in range(1,w):
            j=idx[:,k]; cnt=Dn[:,:s.D].sum(1); pj=np.zeros(M)
            for t in range(s.K):
                dd=s.Hp[j,t]; inn=Dn[ar,dd]; pj+=inn*(s.r[j]/s.rsum[dd])
            pj=pj/np.maximum(cnt,1); bad|=(pj<=0)
            lq+=np.log(np.where(pj>0,pj,1.0)); s._addDn(Dn,j)
        lq[bad]=-np.inf
        return lq
    def logq_pl(s,idx):
        rr=s.r[idx]; cum=np.cumsum(rr,axis=1)-rr
        return (np.log(rr)-np.log(s.R-cum)).sum(1)
    def fail(s,idx):
        M,w=idx.shape; dets=np.zeros((M,s.D+1),dtype=np.uint8); ob=np.zeros(M,dtype=np.uint8); ar=np.arange(M)[:,None]
        for k in range(w):
            j=idx[:,k]; dets[ar,s.Hp[j]]^=1; ob^=s.obs[j]
        pred=s.m.decode_batch(dets[:,:s.D])
        return (pred[:,0]!=ob)
    def batch(s,w,M):
        """returns (sum fail*W, sum (fail*W)^2, M, nfail)"""
        nw=int(s.rng.binomial(M,ALPHA)); parts=[]
        iw=s.sample_walk(nw,w) if nw>0 else np.empty((0,w),dtype=np.int64)
        ip=s.sample_pl(M-nw,w) if M-nw>0 else np.empty((0,w),dtype=np.int64)
        idx=np.concatenate([iw,ip]); 
        srt=np.sort(idx,axis=1); dup=(srt[:,1:]==srt[:,:-1]).any(1) if w>1 else np.zeros(len(idx),bool)
        ok=~dup
        S1=S2=0.0; nf=0
        if ok.any():
            ii=idx[ok]; lt=s.logr[ii].sum(1)-s.loge[w]-gammaln(w+1)
            lw=s.logq_walk(ii); lp=s.logq_pl(ii)
            lmix=np.logaddexp(math.log(ALPHA)+lw, math.log(1-ALPHA)+lp)
            W=np.exp(lt-lmix); f=s.fail(ii)
            v=W[f]; S1=float(v.sum()); S2=float((v**2).sum()); nf=int(f.sum())
        return S1,S2,M,nf
    # ---- main
    def run(s,tag='',rel=0.02,M0=5000,tmax=1500,maxM=40_000_000,log=print):
        d=s.d; w0=(d+1)//2; t0=time.time()
        T=np.cumsum(s.Pw[::-1])[::-1]  # T[w]=P(W>=w)
        Winit=int(np.argmax(s.Pw<1e-40*s.Pw.max()*0+0)) if False else None
        # strata up to where tail P(W>w) < 1e-30
        wl=[w for w in range(1,s.maxw+1)]
        Winit=next((w for w in wl if T[w+1]<1e-30) , s.maxw) if s.maxw+1<len(T) else s.maxw
        st={w:[0.0,0.0,0,0] for w in range(1,Winit+1)}
        for w in st:
            M=M0*(20 if w<w0 else 1); a=s.batch(w,M); st[w]=[a[0],a[1],a[2],a[3]]
        Wmax=1/(1-ALPHA)
        def summary():
            mu=0; U2=0; U={}
            for w,(S1,S2,M,nf) in st.items():
                f=S1/M; mu+=s.Pw[w]*f
                if w<w0 and nf==0: u=0.0
                else:
                    ess=(S1*S1/S2) if S2>0 else 0.0   # effective number of failing samples
                    u=s.Pw[w]*math.sqrt((S2+Wmax**2)/M**2)*max(1.0,math.sqrt(30.0/max(ess,1.0)))
                U[w]=u; U2+=u*u
            return mu,math.sqrt(U2),U
        it=0
        while True:
            it+=1; mu,se,U=summary(); totM=sum(v[2] for v in st.values())
            wcap=max(st); tl=s.tail(wcap)
            if mu>0 and se<=rel*mu and tl<=1e-3*mu: break
            if time.time()-t0>tmax or totM>maxM: log(f'{tag} STOP by budget'); break
            if mu>0 and tl>1e-3*mu:
                pass  # tail covered by Winit<... ; flagged below
            mx=max(U.values()) if U else 0
            sel=[w for w,u in sorted(U.items(),key=lambda x:-x[1]) if u>=0.3*mx and u>0][:6]
            if not sel: break
            for w in sel:
                Mw=st[w][2]; B=int(min(max(M0,Mw),400000))
                a=s.batch(w,B)
                for i in range(4): st[w][i]+=a[i]
            if it%10==0: log(f'{tag} it{it} mu={mu:.4e} relse={se/max(mu,1e-300):.3f} M={totM} t={time.time()-t0:.0f}s sel={sel}')
        mu,se,U=summary(); totM=sum(v[2] for v in st.values()); tl=s.tail(max(st))
        zero_bound=sum(s.Pw[w]*3*Wmax/st[w][2] for w in st if st[w][3]==0 and w>=w0)
        return dict(d=d,p=s.p,truth=mu,se=se,tail=tl,zero_stratum_bound=zero_bound,totM=totM,secs=time.time()-t0,
                    strata={str(w):v for w,v in st.items()},N=s.N,lam=s.lam)
if __name__=='__main__':
    d=int(sys.argv[1]); p=float(sys.argv[2]); seed=int(sys.argv[3]); out=sys.argv[4]
    tmax=float(sys.argv[5]) if len(sys.argv)>5 else 1500
    kind=sys.argv[6] if len(sys.argv)>6 else 'full'
    m=Model(d,p,seed,kind); r=m.run(tag=f'd{d}p{p}',tmax=tmax); r['kind']=kind
    json.dump(r,open(out,'w')); print('RESULT',d,p,r['truth'],r['se'],r['tail'],r['zero_stratum_bound'],r['totM'],round(r['secs']))
