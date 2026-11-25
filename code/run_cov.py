# Frozen-run driver for the L coverage table on a compute host (see PLAN.md "FREEZE ADDENDUM"). Resumable: every finished task is
# appended to OUT/manifest.jsonl keyed by (stage,key); a restart skips keys already present. Run from this directory:
#   PYTHONPATH=../pylib OMP_NUM_THREADS=1 python3 run_cov.py OUT NPROC
# Stages in priority order: S1 main replicate window data (4500 tasks), S2 truthextra (12 real-circuit truth cells, 2nd independent run),
# S3 control A (data-only noise), S4 control B (edge-iid surrogate, 300 reps).
import sys,os,json,time,subprocess
from multiprocessing import Pool
sys.argv_saved=list(sys.argv)
OUT=sys.argv[1]; NPROC=int(sys.argv[2]); os.makedirs(OUT,exist_ok=True)
CH=json.load(open('../results/chunks.json'))
TARGETS=[1e-3,1e-4,1e-5,1e-6]
def t_s1(a):
    import fitdata; return fitdata.one(a)
def t_s2(a):
    d,p,seed=a; f=f'{OUT}/truthextra_d{d}_p{p:g}.json'
    r=subprocess.run(['python3','truth.py',str(d),f'{p:g}',str(seed),f,'2400'],capture_output=True,text=True)
    return dict(d=d,p=p,seed=seed,rc=r.returncode,tail=r.stdout[-300:])
def t_s3(a):
    import control_a; return control_a.one(a)
def t_s4(a):
    import control_b; return control_b.one(a)
FN={'S1':t_s1,'S2':t_s2,'S3':t_s3,'S4':t_s4}
def work(x):
    st,key,a=x; t=time.time(); r=FN[st](a); return dict(stage=st,key=key,res=r,secs=round(time.time()-t,1))
def tasks():
    T=[]
    for r in range(300):
        for d in (3,5,7):
            for ip in range(5): T.append(('S1',f'{r},{d},{ip}',(r,d,ip,CH[f'{d},{ip}'],'full')))
    i=0
    for d in (7,5,3):
        for p in TARGETS: T.append(('S2',f'{d},{p:g}',(d,p,5000+i))); i+=1
    for d in (7,5,3):
        for ip in range(5): T.append(('S3',f'{d},{ip}',(d,ip,3e8)))
    for r in range(300):
        for d in (3,5,7):
            for ip in range(5): T.append(('S4',f'{r},{d},{ip}',(r,d,ip)))
    return T
if __name__=='__main__':
    mf=f'{OUT}/manifest.jsonl'; done=set()
    if os.path.exists(mf):
        for l in open(mf):
            try: x=json.loads(l); done.add((x['stage'],x['key']))
            except Exception: pass
    T=[t for t in tasks() if (t[0],t[1]) not in done]
    print(f'{time.ctime()} total {len(tasks())} done {len(done)} todo {len(T)}',flush=True)
    with Pool(NPROC) as P, open(mf,'a') as f:
        for x in P.imap(work,T,chunksize=1):
            f.write(json.dumps(x)+'\n'); f.flush()
    print(f'{time.ctime()} ALLDONE',flush=True); open(f'{OUT}/ALLDONE','w').write('ok')
