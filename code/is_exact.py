"""IS validation runner (job 4, pre-registered in PLAN.md 'IS validation'). Run from the bundle root (outputs under work/L_is).
Estimator = truth.Model (weight-stratified mixture IS) pointed at the EXACT-enumeration model of exact/common.py.
usage: is_exact.py D P NSEEDS NPROC [TMAX_SEC]   (resumable; per-run json in out/). NEVER used for d=9.
       is_exact.py agg D P   -> per-cell verdict vs work/L_exact/out/dD/manifest.json"""
import sys, os, json, time, math
sys.path.insert(0, os.path.expanduser("work/L_exact")); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, truth
from common import build, lin_dem
OUT = os.path.expanduser("work/L_is/out"); _M = {}
class _Circ:
    def __init__(s, d, p): s.d, s.p = d, p
    def detector_error_model(s, decompose_errors=False):
        M = _M[s.d]
        return M["dem_dec"] if decompose_errors else lin_dem(M, s.p)
def make_model(d, p, seed):
    _M.setdefault(d, build(d)); old = truth.circuit
    truth.circuit = lambda dd, pp, kind="full": _Circ(dd, pp)
    try: return truth.Model(d, p, seed)
    finally: truth.circuit = old
def one(a):
    d, p, i, tmax = a; f = f"{OUT}/d{d}_p{p:g}_s{i:02d}.json"
    if os.path.exists(f): return f
    seed = int(np.random.SeedSequence([4, d, int(round(-math.log10(p))), i]).generate_state(1)[0])
    m = make_model(d, p, seed); r = m.run(tag=f"d{d}p{p:g}s{i}", rel=0.02, tmax=tmax, log=lambda *x: None)
    r.update(seed_index=i, seed=seed); json.dump(r, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); return f
def poly(d, p):
    c = json.load(open(os.path.expanduser(f"work/L_exact/out/d{d}/manifest.json")))["coefficients_p_L=sum c_k p^k"]
    return sum(v["coef"] * p ** int(k) for k, v in c.items())
def agg(d, p):
    ex = poly(d, p); rs = [json.load(open(f"{OUT}/{x}")) for x in sorted(os.listdir(OUT)) if x.startswith(f"d{d}_p{p:g}_s") and x.endswith(".json")]
    ins = [abs(r["truth"] - ex) <= 1.96 * r["se"] for r in rs]; rw = [2 * 1.96 * r["se"] / r["truth"] for r in rs]
    z = [(r["truth"] - ex) / r["se"] for r in rs]
    ok = len(rs) == 20 and sum(ins) >= 18 and float(np.median(rw)) <= 0.20
    print(json.dumps(dict(d=d, p=p, exact=ex, n_runs=len(rs), inside=sum(ins), median_rel_width=float(np.median(rw)), max_rel_width=max(rw), mean_z=float(np.mean(z)), median_secs=float(np.median([r["secs"] for r in rs])), verdict="PASS" if ok else "FAIL (or incomplete)")))
if __name__ == "__main__":
    if sys.argv[1] == "agg": agg(int(sys.argv[2]), float(sys.argv[3])); sys.exit()
    from multiprocessing import Pool
    d = int(sys.argv[1]); p = float(sys.argv[2]); n = int(sys.argv[3]); NP = int(sys.argv[4]); tmax = float(sys.argv[5]) if len(sys.argv) > 5 else 1500
    assert d <= 5, "d>5 refused: holdout must not run here"
    os.makedirs(OUT, exist_ok=True)
    with Pool(NP) as pool:
        for f in pool.imap_unordered(one, [(d, p, i, tmax) for i in range(n)]): print("done", f, flush=True)
    agg(d, p)
