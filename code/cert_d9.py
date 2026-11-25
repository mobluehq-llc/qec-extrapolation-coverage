"""d=9, p=1e-4 CERTIFIED interval: a NEW, separately pre-registered measurement (NOT a rerun of the one-shot holdout; PLAN.md
'Certified interval amendment'). Sampling allocation is fixed from the OLD holdout strata data (independent of stage 2);
only fresh stage-2 samples (new seeds) enter the bounds.
  alloc                       -> ../results/cert_d9_alloc.json (committed before the run)
  work MACHINE NMACH NPROC OUT -> process tasks i with i%NMACH==MACHINE (resumable)
  combine OUT                 -> ../results/cert_d9_result.json"""
import sys, os, json, glob, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import certified as C
D, P = 9, 1e-4
NTOT = 10**9; NFREE = 10**7; CH = 5_000_000; A0 = {1, 2, 3, 4}
def holdout():
    rs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(HERE, "..", "results", "holdout_truth_d9_p1e-4_s*.json")))]
    agg = {}
    for r in rs:
        for w, (S1, S2, M, nf) in r["strata"].items():
            a = agg.setdefault(int(w), [0.0, 0.0, 0]); a[0] += S1; a[1] += S2; a[2] += M
    return agg, rs
def plan():
    agg, rs = holdout(); m = C.make_model("full", D, P, 1)
    phat = sum(m.Pw[w] * agg[w][0] / agg[w][2] for w in agg if w >= 5 and w <= 40)
    wmax = next(w for w in range(5, 60) if C.tail_mass(m, w) <= 1e-3 * phat)
    strata = list(range(5, wmax + 1)); alloc = C.allocate(m, strata, set(), {w: tuple(agg[w]) for w in strata}, NTOT)
    for w in (1, 2, 3, 4): alloc[w] = NFREE
    return m, dict(d=D, p=P, A0=sorted(A0), Ntot=NTOT, nfree=NFREE, chunk=CH, wmax=wmax, phat_holdout_prior=phat, tail=C.tail_mass(m, wmax),
        alloc={str(w): int(n) for w, n in sorted(alloc.items())}, b={str(w): C.bw(m, w) for w in range(1, wmax + 1)},
        P_w={str(w): float(m.Pw[w]) for w in range(1, wmax + 1)}, N=m.N, dem_lam=float(m.lam))
def tasks():
    pl = json.load(open(os.path.join(HERE, "..", "results", "cert_d9_alloc.json"))); T = []
    for w, n in pl["alloc"].items():
        w = int(w); k = 0
        while n > 0: c = min(CH, n); T.append((w, k, c)); n -= c; k += 1
    return pl, T
def one(a):
    w, k, c, out = a; f = f"{out}/w{w}_c{k}.npz"
    if os.path.exists(f): return f
    seed = int(np.random.SeedSequence([8, D, w, k]).generate_state(1)[0]); m = C.make_model("full", D, P, seed)
    v, n = C.draw(m, w, c); np.savez(f + ".tmp.npz", v=v, n=n); os.replace(f + ".tmp.npz", f); return f
def combine(out):
    pl, T = tasks(); m = C.make_model("full", D, P, 1); data = {}
    for (w, k, c) in T:
        z = np.load(f"{out}/w{w}_c{k}.npz"); assert int(z["n"]) == c
        v, n = data.get(w, (np.empty(0), 0)); data[w] = (np.concatenate([v, z["v"]]), n + c)
    wmax = pl["wmax"]; strata = list(range(1, wmax + 1))
    a0 = C.combine(m, strata, data, set(A0)); free = C.combine(m, strata, data, set())
    res = dict(plan=pl, A0=a0, assumption_free=free, d=D, p=P, nsamples={w: int(data[w][1]) for w in data}, nfail={w: int(len(data[w][0])) for w in data})
    json.dump(res, open(os.path.join(HERE, "..", "results", "cert_d9_result.json"), "w"), indent=1, default=float)
    print("A-0 interval", a0["L"], a0["U"], "est", a0["est"]); print("assumption-free", free["L"], free["U"])
if __name__ == "__main__":
    c = sys.argv[1]
    if c == "alloc":
        m, pl = plan(); json.dump(pl, open(os.path.join(HERE, "..", "results", "cert_d9_alloc.json"), "w"), indent=1); print(json.dumps(pl, indent=1))
    elif c == "work":
        from multiprocessing import Pool
        mach, nm, npr, out = {int(x) for x in sys.argv[2].split(',')}, int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; os.makedirs(out, exist_ok=True)
        pl, T = tasks(); mine = [(w, k, c_, out) for i, (w, k, c_) in enumerate(T) if i % nm in mach]; print(len(mine), "tasks", flush=True)
        with Pool(npr) as pool:
            for f in pool.imap_unordered(one, mine): print("done", f, flush=True)
    elif c == "combine": combine(sys.argv[2])
