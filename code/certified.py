"""CERTIFIED p_L interval (see PLAN.md 'Certified interval amendment').

p_L = sum_w P(W=w) f_w,  f_w in [0,1],  P(W=w) exact (elementary symmetric polynomials of the odds).
Per stratum w the weight-stratified mixture IS (truth.Model samplers, UNCHANGED) gives iid Y = fail*Wt in [0,b_w],
  Wt = t/mix <= t/((1-ALPHA) q_PL) <= B_w/(1-ALPHA),  B_w = R^w/(e_w w!)   (R = sum odds; R-cum <= R)
  so b_w = B_w/(1-ALPHA) is a RIGOROUS range bound (not the heuristic 1/(1-ALPHA) quoted in truth.py).
Each stratum mean f_w gets a two-sided betting confidence interval (Waudby-Smith & Ramdas 2023 style, fixed lambda grid
  averaged => test martingale; Ville/Markov), at level delta/(2 S) each (union over 2S one-sided events).
  L = sum P_w LCB_w ;  U = sum P_w UCB_w + P(W>Wmax)  (f<=1 in the tail).  Strata known to have f_w=0 by COMPLETE
  exhaustive enumeration (exact/ manifests) or by a labelled assumption (A-0) are fixed at 0.
Stage design: pilot -> allocation fixed -> FRESH stage-2 samples only enter the bounds (pilot discarded)."""
import sys, os, json, math, time
import numpy as np
from scipy.special import gammaln, logsumexp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import truth
ALPHA = truth.ALPHA
DELTA = 0.05

def make_model(kind, d, p, seed):
    if kind == "exact":
        import is_exact
        return is_exact.make_model(d, p, seed)
    return truth.Model(d, p, seed, "full")

def bw(m, w):
    """rigorous range bound for Y=fail*Wt in stratum w"""
    return math.exp(w * math.log(m.R) - m.loge[w] - gammaln(w + 1)) / (1 - ALPHA)

def sample_vals(m, w, M):
    """M iid mixture draws in stratum w; returns (nonzero values of fail*Wt, max Wt over all valid draws)."""
    nw = int(m.rng.binomial(M, ALPHA))
    iw = m.sample_walk(nw, w) if nw > 0 else np.empty((0, w), dtype=np.int64)
    ip = m.sample_pl(M - nw, w) if M - nw > 0 else np.empty((0, w), dtype=np.int64)
    idx = np.concatenate([iw, ip])
    srt = np.sort(idx, axis=1); dup = (srt[:, 1:] == srt[:, :-1]).any(1) if w > 1 else np.zeros(len(idx), bool)
    ok = ~dup; wmax = 0.0; v = np.empty(0)
    if ok.any():
        ii = idx[ok]; lt = m.logr[ii].sum(1) - m.loge[w] - gammaln(w + 1)
        lw = m.logq_walk(ii); lp = m.logq_pl(ii)
        lmix = np.logaddexp(math.log(ALPHA) + lw, math.log(1 - ALPHA) + lp)
        W = np.exp(lt - lmix); f = m.fail(ii); v = W[f]; wmax = float(W.max())
    return v, wmax

# ---------------- betting confidence bounds for the mean of Y in [0,b] (n draws, nonzero values v, rest are 0)
GP = np.geomspace(1e-3, 30.0, 18)     # lambda/(1/b) multipliers for the 'mean > m' test (K+)
GM = np.geomspace(1e-3, 0.999, 18)    # for the 'mean < m' test (K-), needs lambda (b-m) < 1

def _logK(v, n, b, m0, plus):
    n0 = n - len(v); g = GP if plus else GM; lam = g / b; out = np.full(len(lam), -np.inf)
    for j, l in enumerate(lam):
        if plus:
            if l * m0 >= 1 - 1e-12: continue
            s = n0 * math.log1p(-l * m0) + (np.log1p(l * (v - m0)).sum() if len(v) else 0.0)
        else:
            if l * (b - m0) >= 1 - 1e-12: continue
            s = n0 * math.log1p(l * m0) + (np.log1p(l * (m0 - v)).sum() if len(v) else 0.0)
        out[j] = s
    return float(logsumexp(out) - math.log(len(lam)))

def bounds(v, n, b, dprime):
    """(LCB, UCB) each valid with prob >= 1-dprime. Monotone bisection on log m."""
    thr = math.log(1.0 / dprime); v = np.asarray(v, dtype=float)
    if len(v) and v.max() > b * (1 + 1e-9): raise AssertionError(f"range bound violated: {v.max()} > {b}")
    # LCB = inf{m : K+(m) < 1/dprime}; K+ nonincreasing in m
    lo, hi = -30.0, math.log10(b)
    if _logK(v, n, b, 10 ** lo, True) < thr: lcb = 0.0
    else:
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            if _logK(v, n, b, 10 ** mid, True) < thr: hi = mid
            else: lo = mid
        lcb = 10 ** hi
    # UCB = sup{m : K-(m) < 1/dprime}; K- nondecreasing in m
    lo, hi = -30.0, math.log10(b)
    if _logK(v, n, b, b * (1 - 1e-9), False) < thr: ucb = b
    else:
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            if _logK(v, n, b, 10 ** mid, False) < thr: lo = mid
            else: hi = mid
        ucb = 10 ** hi
    return lcb, ucb

def exact_zero_strata(d):
    """weights whose shell was enumerated COMPLETELY with 0 failing sets (exact/ manifests) -> f_w = 0 rigorously"""
    f = os.path.join(HERE, "..", "exact", f"d{d}_manifest.json")
    if not os.path.exists(f): f = os.path.expanduser(f"work/L_exact/out/d{d}/manifest.json")
    if not os.path.exists(f): return set()
    sh = json.load(open(f))["shells"]
    return {int(w) for w, s in sh.items() if s.get("complete") and s.get("n_fail") == 0}

def strata_range(m):
    T = np.cumsum(m.Pw[::-1])[::-1]
    return next((w for w in range(1, m.maxw + 1) if T[w + 1] < 1e-30), m.maxw)

def tail_mass(m, wmax): return float(m.Pw[wmax + 1:].sum())

def combine(m, strata, data, zero, delta=DELTA, extra_upper=0.0):
    """data[w] = (v, n). zero = strata fixed at f=0. returns dict with L,U,est"""
    act = [w for w in strata if w not in zero]; dp = delta / (2 * len(act)); L = U = est = 0.0; per = {}
    for w in act:
        v, n = data[w]; b = bw(m, w); l, u = bounds(v, n, b, dp)
        L += m.Pw[w] * l; U += m.Pw[w] * u; est += m.Pw[w] * (v.sum() / n)
        per[w] = dict(P=float(m.Pw[w]), n=int(n), nfail=int(len(v)), f_hat=float(v.sum() / n), lcb=l, ucb=u, b=b)
    wmax = max(strata); U += tail_mass(m, wmax) + extra_upper
    return dict(L=L, U=U, est=est, tail=tail_mass(m, wmax), per=per, n_events=2 * len(act), delta=delta)

def allocate(m, strata, zero, S1S2M, Ntot, nmin=2000):
    """S1S2M[w]=(S1,S2,M) from a PILOT / prior data (independent of stage 2). n_w ~ P_w sqrt((S2+b^2)/M)"""
    a = {}
    for w in strata:
        if w in zero: continue
        S1, S2, M = S1S2M[w]; b = bw(m, w); a[w] = m.Pw[w] * math.sqrt((S2 + b * b) / M)
    tot = sum(a.values()); return {w: max(nmin, int(round(Ntot * a[w] / tot))) for w in a}

def draw(m, w, n, chunk=20000):
    vs = []; done = 0
    while done < n:
        k = min(chunk, n - done); v, _ = sample_vals(m, w, k); vs.append(v); done += k
    return (np.concatenate(vs) if vs else np.empty(0)), n

def replicate(kind, d, p, rep, Ntot, n0=20000):
    """full self-contained certified interval for one replicate (d=3/5 validation)."""
    seed = int(np.random.SeedSequence([7, d, int(round(-math.log10(p))), rep]).generate_state(1)[0])
    m = make_model(kind, d, p, seed); t0 = time.time(); zero = exact_zero_strata(d) if kind == "exact" else set()
    wtop = strata_range(m); pil = {}
    for w in range(1, wtop + 1):
        if w in zero: continue
        v, n = draw(m, w, n0); pil[w] = (float(v.sum()), float((v ** 2).sum()), n)
    phat = sum(m.Pw[w] * pil[w][0] / pil[w][2] for w in pil)
    wmax = next((w for w in range(1, wtop + 1) if tail_mass(m, w) <= 1e-3 * max(phat, 1e-300)), wtop)
    strata = list(range(1, wmax + 1)); alloc = allocate(m, strata, zero, pil, Ntot)
    data = {w: draw(m, w, alloc[w]) for w in alloc}
    r = combine(m, strata, data, zero); r.update(kind=kind, d=d, p=p, rep=rep, seed=seed, wmax=wmax, zero=sorted(zero), Ntot=Ntot,
        alloc={str(w): n for w, n in alloc.items()}, secs=time.time() - t0, phat_pilot=phat); return r

def truth_poly(d, p):
    f = os.path.join(HERE, "..", "exact", f"d{d}_manifest.json")
    if not os.path.exists(f): f = os.path.expanduser(f"work/L_exact/out/d{d}/manifest.json")
    c = json.load(open(f))["coefficients_p_L=sum c_k p^k"]; return sum(v["coef"] * p ** int(k) for k, v in c.items())

def _job(a):
    kind, d, p, rep, Ntot, out = a; f = f"{out}/cert_{kind}_d{d}_p{p:g}_r{rep:02d}.json"
    if os.path.exists(f): return f
    r = replicate(kind, d, p, rep, Ntot); json.dump(r, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); return f

def agg(kind, d, p, out, nmin=20):
    ex = truth_poly(d, p); rs = [json.load(open(f"{out}/{x}")) for x in sorted(os.listdir(out)) if x.startswith(f"cert_{kind}_d{d}_p{p:g}_r")]
    ins = [r["L"] <= ex <= r["U"] for r in rs]; rw = [(r["U"] - r["L"]) / ex for r in rs]
    ok = len(rs) >= nmin and all(ins) and float(np.median(rw)) <= 0.50
    print(json.dumps(dict(d=d, p=p, exact=ex, n=len(rs), inside=int(sum(ins)), median_rel_width=float(np.median(rw)), max_rel_width=max(rw),
        median_L_over_truth=float(np.median([r["L"] / ex for r in rs])), median_U_over_truth=float(np.median([r["U"] / ex for r in rs])),
        median_secs=float(np.median([r["secs"] for r in rs])), coverage_ok=bool(len(rs) >= nmin and all(ins)), width_ok=bool(float(np.median(rw)) <= 0.5),
        verdict="PASS" if ok else "FAIL")))

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "agg": agg(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]); sys.exit()
    if cmd == "rep":   # rep KIND D P NREPS NPROC NTOT OUTDIR
        from multiprocessing import Pool
        kind, d, p, nrep, NP, Ntot, out = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(float(sys.argv[7])), sys.argv[8]
        assert kind == "exact" and d <= 5, "validation runs only (exact model, d<=5)"
        os.makedirs(out, exist_ok=True)
        with Pool(NP) as pool:
            for f in pool.imap_unordered(_job, [(kind, d, p, i, Ntot, out) for i in range(nrep)]): print("done", f, flush=True)
        agg(kind, d, p, out)
