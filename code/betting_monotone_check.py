"""Check of the lower (and upper) betting inversion in certified.py on the REAL d=9 stage-2 data (results/cert_d9_chunks).
Question (ChatGPT final check on v7): the admissible lambda grid depends on the trial mean m0 (lambda*m0 < 1 for the K+ test);
is the root search monotone? Facts checked: (1) K+(m) and K-(m) evaluated on a dense log-m grid are non-increasing / non-decreasing
(tolerance 1e-9 in log K); (2) the first cutoff m = b/g_max and the number of admissible multipliers at the computed bounds;
(3) brute-force root on the dense grid equals the bisection LCB/UCB of certified.bounds to the grid spacing. Run on a compute host:  python3 code/betting_monotone_check.py"""
import os, sys, json, math, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import certified as C
R = os.path.join(HERE, "..", "results"); pl = json.load(open(f"{R}/cert_d9_alloc.json")); D = 16  # 16 one-sided statements under A-0 -> delta' = 0.05/16
dprime = C.DELTA / D; thr = math.log(1/dprime); out = {}; worst_mono = 0.0; worst_jump = 0.0
for w in range(5, 13):
    v = np.concatenate([np.load(f)["v"] for f in sorted(glob.glob(f"{R}/cert_d9_chunks/w{w}_c*.npz"))]); n = int(sum(int(np.load(f)["n"]) for f in glob.glob(f"{R}/cert_d9_chunks/w{w}_c*.npz"))); b = pl["b"][str(w)]
    lcb, ucb = C.bounds(v, n, b, dprime); fhat = v.sum()/n
    ms = np.geomspace(max(fhat*1e-4, 1e-30), b*0.999, 4000)
    kp = np.array([C._logK(v, n, b, m, True) for m in ms]); km = np.array([C._logK(v, n, b, m, False) for m in ms])
    fin = np.isfinite(kp); dkp = np.diff(kp[fin]); dkm = np.diff(km[np.isfinite(km)])
    mono = max(float(dkp.max(initial=-1e9)), -float(dkm.min(initial=1e9)))   # >0 would be a violation
    # smallest m at which any multiplier becomes inadmissible (lambda_j m >= 1, i.e. m >= b/g_max) and the admissible count AT the bounds
    cut = b / float(C.GP.max()); nadm = lambda m: int(sum(1 for g in C.GP if (g/b)*m < 1 - 1e-12)); adm_roots = min(nadm(max(lcb, 1e-300)), nadm(ucb))
    # admissible-grid count along m
    adm = [int(sum(1 for g in C.GP if (g/b)*m < 1 - 1e-12)) for m in ms[::400]]
    # brute-force roots
    below = np.where(kp < thr)[0]; bf_l = float(ms[below[0]]) if len(below) else 0.0
    aboveU = np.where(km >= thr)[0]; bf_u = float(ms[aboveU[0]]) if len(aboveU) else b
    spacing = float(np.exp(np.log(ms[1]/ms[0])))-1
    out[w] = dict(n=n, nfail=int(len(v)), b=b, lcb=lcb, ucb=ucb, bf_lcb=bf_l, bf_ucb=bf_u, rel_diff_lcb=(abs(bf_l-lcb)/lcb if lcb > 0 else None), rel_diff_ucb=abs(bf_u-ucb)/ucb, grid_step=spacing,
                  logKplus_nonincreasing_violation=float(max(0.0, float(dkp.max(initial=-1e9)))), logKminus_nondecreasing_violation=float(max(0.0, -float(dkm.min(initial=1e9)))), first_cutoff_m=cut, cutoff_over_max_bound=cut/max(lcb, ucb, 1e-300), admissible_at_bounds=adm_roots, admissible_counts_along_m=adm)
    worst_mono = max(worst_mono, out[w]["logKplus_nonincreasing_violation"], out[w]["logKminus_nondecreasing_violation"]); worst_jump = min(worst_jump or 1e9, out[w]['admissible_at_bounds'])
    print(w, {k: (round(x, 4) if isinstance(x, float) else x) for k, x in out[w].items() if k != "admissible_counts_along_m"}, out[w]["admissible_counts_along_m"], flush=True)
json.dump(dict(strata=out, worst_monotonicity_violation=worst_mono, min_admissible_at_bounds=worst_jump, delta_prime=dprime), open(f"{R}/L_betting_monotone_v7.json", "w"), indent=1)
print("worst monotonicity violation", worst_mono, "min admissible multipliers at any bound (of 18)", worst_jump)
