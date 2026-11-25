# Below-threshold amendment driver (resumable). Run from this directory: PYTHONPATH=../pylib OMP_NUM_THREADS=1 python3 run_bw.py OUT NPROC
import sys,os,json,time
from multiprocessing import Pool
import fitdata_bw
OUT=sys.argv[1]; NPROC=int(sys.argv[2]); os.makedirs(OUT,exist_ok=True)
CH=json.load(open('../results/chunks_bw.json'))
def work(a):
    t=time.time(); r=fitdata_bw.one(a); return dict(key=f'{a[0]},{a[1]},{a[2]}',res=r,secs=round(time.time()-t,1))
if __name__=='__main__':
    mf=f'{OUT}/manifest_bw.jsonl'; done=set()
    if os.path.exists(mf):
        for l in open(mf):
            try: done.add(json.loads(l)['key'])
            except Exception: pass
    T=[(r,d,ip,CH[f'{d},{ip}'],'full') for ip in range(5) for d in (7,5,3) for r in range(300)]   # costly cells first
    T=[a for a in T if f'{a[0]},{a[1]},{a[2]}' not in done]
    print(time.ctime(),'todo',len(T),flush=True)
    with Pool(NPROC) as P, open(mf,'a') as f:
        for x in P.imap_unordered(work,T,chunksize=1): f.write(json.dumps(x)+'\n'); f.flush()
    rows=[json.loads(l)['res'] for l in open(mf)]; json.dump(sorted(rows,key=lambda x:(x['rep'],x['d'],x['ip'])),open(f'{OUT}/fit_bw.json','w')); print(time.ctime(),'ALLDONE',len(rows),flush=True)
