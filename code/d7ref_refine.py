"""L d=7 reference refinement (PLAN.md, Gemini amendment item 3; pre-registered, no adaptation).
Run from this directory with PYTHONPATH=<pylib> OMP_NUM_THREADS=1.
  alloc            -> ../results/L_d7ref_alloc.json  (pooled pilot = the 8 stored d=7 runs, which are DISCARDED from the estimate)
  work NPROC OUT   -> resumable fixed-allocation sampling; chunk files OUT/p{p}_w{w}_c{k}.json hold S1,S2,M,nf (fresh seeds, SeedSequence prefix 16)
  combine OUT      -> ../results/L_d7ref_result.json : per cell the SUM of per-stratum means (pooled S1/M), plain se and frozen conservative se.
No stopping on data, no 1/sigma^2 weights. Strata w<=3: 1e5 draws, treated as 0 (A-0 type), as in the frozen runs."""
import sys, os, json, glob, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "..", "results"); sys.path.insert(0, HERE)
PS = ["1e-3", "1e-4", "1e-5", "1e-6"]; TARGET = 0.015; SAFETY = 1.5; FLOOR = 100000; CH = 1_000_000; SUB = 100_000; W0 = 4; ALPHA = 0.7
D = 7
def pilot(ps):
    fs = [os.path.join(RES, f"truth_d7_p{ps}.json")] + [f for f in glob.glob(os.path.join(RES, "truthextra_d7_p*.json"))
          if abs(float(os.path.basename(f)[len("truthextra_d7_p"):-5]) - float(ps)) < 1e-9 * float(ps)]
    agg = {}
    for f in fs:
        for w, (S1, S2, M, nf) in json.load(open(f))["strata"].items():
            a = agg.setdefault(int(w), [0.0, 0.0, 0, 0]); a[0] += S1; a[1] += S2; a[2] += M; a[3] += nf
    return agg, len(fs)
def alloc():
    import truth
    out = dict(design=f"Neyman M_w ~ P_w sigma_w from pooled pilot, x{SAFETY} safety, floor {FLOOR}, target plain rel se {TARGET}; w<=3 fixed {FLOOR} draws treated 0", chunk=CH, sub=SUB, cells={})
    for ps in PS:
        agg, nf = pilot(ps); m = truth.Model(D, float(ps), 1, "full")
        P = {w: float(m.Pw[w]) for w in agg}; sig = {}
        for w, (S1, S2, M, _) in agg.items():
            sig[w] = math.sqrt(max(S2 / M - (S1 / M) ** 2, 0.0)) if M > 0 and w >= W0 else 0.0
        mu = sum(P[w] * agg[w][0] / agg[w][2] for w in agg if w >= W0); tot = sum(P[w] * sig[w] for w in agg)
        K = tot / (TARGET * mu) ** 2
        M_w = {w: (FLOOR if w < W0 else int(max(FLOOR, math.ceil(SAFETY * K * P[w] * sig[w])))) for w in sorted(agg)}
        out["cells"][ps] = dict(pilot_runs=nf, pilot_mu=mu, K=K, sum_Psig=tot, P_w={str(w): P[w] for w in P}, sigma_w={str(w): sig[w] for w in sig},
                                M_w={str(w): M_w[w] for w in M_w}, total_draws=sum(M_w.values()))
        print(ps, "pilot runs", nf, "mu", mu, "total draws %.3e" % sum(M_w.values()), flush=True)
    out["grand_total_draws"] = sum(c["total_draws"] for c in out["cells"].values())
    json.dump(out, open(os.path.join(RES, "L_d7ref_alloc.json"), "w"), indent=1); print("grand total %.3e" % out["grand_total_draws"])
def tasks():
    pl = json.load(open(os.path.join(RES, "L_d7ref_alloc.json"))); T = []
    for ci, ps in enumerate(PS):
        for w, n in pl["cells"][ps]["M_w"].items():
            w = int(w); k = 0
            while n > 0: c = min(CH, n); T.append((ci, ps, w, k, c)); n -= c; k += 1
    T.sort(key=lambda t: (t[3], t[0], -t[2]))   # chunk index outermost: balanced partial progress
    return T
_models = {}
def one(a):
    ci, ps, w, k, c, out = a; f = f"{out}/p{ps}_w{w}_c{k}.json"
    if os.path.exists(f): return f
    import truth
    seed = int(np.random.SeedSequence([16, ci, w, k]).generate_state(1)[0])
    if ps not in _models: _models[ps] = truth.Model(D, float(ps), seed, "full")
    m = _models[ps]; m.rng = np.random.default_rng(seed)
    S1 = S2 = 0.0; M = nf = 0; left = c
    while left > 0:
        b = min(SUB, left); a_ = m.batch(w, b); S1 += a_[0]; S2 += a_[1]; M += a_[2]; nf += a_[3]; left -= b
    assert M == c
    json.dump(dict(p=ps, w=w, k=k, S1=S1, S2=S2, M=M, nf=nf, seed=seed), open(f + ".tmp", "w")); os.replace(f + ".tmp", f); return f
def work(nproc, out):
    from multiprocessing import Pool
    os.makedirs(out, exist_ok=True); T = [t + (out,) for t in tasks()]; todo = [t for t in T if not os.path.exists(f"{out}/p{t[1]}_w{t[2]}_c{t[3]}.json")]
    print(time.ctime(), "tasks", len(T), "todo", len(todo), flush=True); n = 0
    with Pool(nproc) as pool:
        for f in pool.imap_unordered(one, todo, chunksize=1):
            n += 1
            if n % 20 == 0: print(time.ctime(), "done", n, "/", len(todo), flush=True)
    print(time.ctime(), "ALLDONE", flush=True)
def combine(out):
    pl = json.load(open(os.path.join(RES, "L_d7ref_alloc.json"))); res = {}
    for ps in PS:
        c = pl["cells"][ps]; st = {}
        for f in glob.glob(f"{out}/p{ps}_w*_c*.json"):
            r = json.load(open(f)); a = st.setdefault(r["w"], [0.0, 0.0, 0, 0]); a[0] += r["S1"]; a[1] += r["S2"]; a[2] += r["M"]; a[3] += r["nf"]
        mu = 0.0; v = 0.0; U2 = 0.0; Wmax = 1 / (1 - ALPHA); missing = []
        for w, n in c["M_w"].items():
            w = int(w); P = c["P_w"][str(w)]
            if w not in st or st[w][2] != n: missing.append(w); continue
            S1, S2, M, nf = st[w]
            if w < W0: continue   # A-0 type: fixed at 0
            f = S1 / M; mu += P * f; v += P * P * max(S2 / M - f * f, 0.0) / M
            ess = (S1 * S1 / S2) if S2 > 0 else 0.0
            u = P * math.sqrt((S2 + Wmax ** 2) / M ** 2) * max(1.0, math.sqrt(30.0 / max(ess, 1.0))) if nf > 0 else 0.0
            U2 += u * u
        res[ps] = dict(estimate=mu, se_plain=math.sqrt(v), rel_se_plain=math.sqrt(v) / mu if mu else None, se_conservative=math.sqrt(U2),
                       rel_se_conservative=math.sqrt(U2) / mu if mu else None, strata_missing_or_incomplete=missing, nfail={str(w): st[w][3] for w in sorted(st)},
                       draws=sum(x[2] for x in st.values()))
        print(ps, res[ps]["estimate"], "rel se plain", res[ps]["rel_se_plain"], "cons", res[ps]["rel_se_conservative"], "missing", len(missing))
    json.dump(res, open(os.path.join(RES, "L_d7ref_result.json"), "w"), indent=1)
if __name__ == "__main__":
    c = sys.argv[1]
    if c == "alloc": alloc()
    elif c == "work": work(int(sys.argv[2]), sys.argv[3])
    elif c == "combine": combine(sys.argv[2])
