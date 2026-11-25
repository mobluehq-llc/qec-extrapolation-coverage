"""A0-STRUCT: distances + weighted structural test for A-0. usage: a0_struct.py D"""
import sys, time, heapq, math, numpy as np, stim, pymatching
D=int(sys.argv[1]); p=1e-4
c=stim.Circuit.generated("surface_code:rotated_memory_z",distance=D,rounds=D,
  after_clifford_depolarization=p,before_measure_flip_probability=p,after_reset_flip_probability=p)
dem=c.detector_error_model(decompose_errors=True)
t=time.time(); g=c.shortest_graphlike_error(); print("d",D,"shortest_graphlike_error len",len(g),round(time.time()-t,1),"s",flush=True)
dem0=c.detector_error_model()
t=time.time()
SEARCH=len(sys.argv)<=2 or sys.argv[2] not in('nosearch',)
try:
    if not SEARCH: raise RuntimeError('search skipped by flag')
    u=c.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4,
        dont_explore_edges_with_degree_above=4,dont_explore_edges_increasing_symptom_degree=False,canonicalize_circuit_errors=True)
    print("search_for_undetectable_logical_errors len",len(u),round(time.time()-t,1),"s",flush=True)
except Exception as e: print("search err",str(e)[:200],flush=True)
if len(sys.argv)>2 and sys.argv[2]=='dist': sys.exit()
m=pymatching.Matching.from_detector_error_model(dem)
ed=m.edges(); nd=dem.num_detectors; B=nd  # boundary vertex
E={}
for u,v,a in ed:
    v=B if v is None or v<0 else v; u=B if u is None or u<0 else u
    E[(min(u,v),max(u,v))]=(a["weight"],1 if 0 in a["fault_ids"] else 0)
ws=[x[0] for x in E.values()]; print("graph edges",len(E),"w min/max",min(ws),max(ws),"ratio",max(ws)/min(ws))
# odd-parity min even subgraph via double cover
adj=[[] for _ in range(nd+1)]
for (u,v),(w,o) in E.items(): adj[u].append((v,w,o)); adj[v].append((u,w,o))
def dij(s):
    dist={(s,0):0.0}; h=[(0.0,s,0)]
    while h:
        d0,v,par=heapq.heappop(h)
        if d0>dist.get((v,par),1e300): continue
        if v==s and par==1: return d0
        for x,w,o in adj[v]:
            k=(x,par^o); nd_=d0+w
            if nd_<dist.get(k,1e300): dist[k]=nd_; heapq.heappush(h,(nd_,x,par^o))
    return 1e300
Wmin=min(dij(s) for s in range(nd+1)); print("Wmin",Wmin,"= edges-equiv",Wmin/min(ws))
# mechanism weights W_i = sum of component graph weights
Wi=[]
for ins in dem.flattened():
    if ins.type!="error": continue
    comp=[[]]
    for t_ in ins.targets_copy():
        if t_.is_separator(): comp.append([])
        elif t_.is_relative_detector_id(): comp[-1].append(t_.val)
    s=0.0
    for cc in comp:
        if len(cc)==0: continue
        k=(min(cc[0],B),max(cc[0],B)) if len(cc)==1 else (min(cc),max(cc))
        if len(cc)>2: s=1e300; continue
        s+=E[k][0]
    Wi.append(s)
Wi=np.sort(np.array(Wi))[::-1]; print("mechanisms",len(Wi),"top W_i",Wi[:5],"max comps")
S4=Wi[:4].sum(); marg=64*max(ws)/2**14
print("S4",S4,"2*S4",2*S4,"margin",marg,"CHEAP TEST (Wmin>2S4+margin):",Wmin>2*S4+marg)
T=Wmin/2
# exact count of subsets (size 1..4) with S>=T : brute count via sorted array + searchsorted
asc=np.sort(Wi); n=len(asc)
cnt={1:int((asc>=T).sum())}
import itertools
def count_ge(k):
    if k==2:
        tot=0
        for i in range(n):
            j=np.searchsorted(asc,T-asc[i],side="left"); j=max(j,i+1); tot+=n-j
        return tot
    if k==3:
        tot=0
        for i in range(n):
            for j in range(i+1,n):
                l=np.searchsorted(asc,T-asc[i]-asc[j],side="left"); l=max(l,j+1); tot+=n-l
        return tot
cnt[2]=count_ge(2)
if n<=2500 or True:
    cnt[3]=count_ge(3) if n<=14000 else None
print("residual subsets with S>=Wmin/2 (mechanism-weight upper bound), exact:",cnt,"total C(n,w)",[math.comb(n,w) for w in (1,2,3)])
