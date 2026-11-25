"""Numbers for draft v3 (method section for the reference estimator; certified-interval check on the surrogate).
Run from this folder after make_numbers.py and make_numbers_amend.py:
  python3 make_numbers_v3.py
Reads only stored result files (../results/truth_*.json, truthextra_*.json, holdout_truth_*.json, cert_validation/, cert_d9_alloc.json).
Nothing is simulated. Writes numbers_v3.tex (keys m.*, c.*) and numbers_v3.json (key -> value, source, how)."""
import json, glob, math, re, os
import numpy as np
R = "../results"; V = {}
def put(k, v, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), src, how)
TG = [1e-3, 1e-4, 1e-5, 1e-6]
runs = {}
for f in sorted(glob.glob(f"{R}/truth_d*_p*.json")) + sorted(glob.glob(f"{R}/truthextra_d*_p*.json")):
    x = json.load(open(f))
    if x["d"] in (3, 5, 7) and any(abs(x["p"] / t - 1) < 1e-9 for t in TG): runs.setdefault(x["d"], []).append(x)
assert {d: len(v) for d, v in runs.items()} == {3: 8, 5: 8, 7: 8}
SRC = "results/truth_d*_p*.json + truthextra_d*_p*.json (24 runs: 12 cells x 2)"
for d in (3, 5, 7):
    v = runs[d]; w0 = (d + 1) // 2
    assert len({x["N"] for x in v}) == 1 and all(x["kind"] == "full" for x in v)
    put(f"m.N.d{d}", v[0]["N"], SRC, "fault mechanisms of the circuit detector error model (field N)")
    assert all(sorted(int(w) for w in x["strata"]) == list(range(1, 121)) for x in v)
    for x in v:
        st = x["strata"]
        assert all(st[str(w)][3] == 0 and st[str(w)][2] == 100000 for w in range(1, w0)), "low strata: 100000 draws, no failure"
        assert st[str(w0)][3] > 0 and x["tail"] == 0.0
    tm = [x["totM"] / 1e6 for x in v]
    put(f"m.totM.d{d}.min", f"{min(tm):.2f}" if d == 3 else f"{min(tm):.1f}", SRC, "samples per run, millions (totM)")
    put(f"m.totM.d{d}.max", f"{max(tm):.2f}" if d == 3 else f"{max(tm):.1f}", SRC, "samples per run, millions (totM)")
    rs = [100 * x["se"] / x["truth"] for x in v]
    put(f"m.runrelse.d{d}.min", f"{min(rs):.1f}", SRC, "single-run se/estimate, percent"); put(f"m.runrelse.d{d}.max", f"{max(rs):.1f}", SRC, "")
    sc = [x["secs"] for x in v]; put(f"m.secs.d{d}.max", f"{max(sc):.0f}", SRC, "run time, s"); put(f"m.secs.d{d}.min", f"{min(sc):.0f}", SRC, "")
low = runs[3] + runs[5]
assert all(x["se"] <= 0.02 * x["truth"] and x["secs"] < 2400 for x in low), "d=3,5 runs stopped by the 2% rule"
assert all(x["se"] > 0.02 * x["truth"] and x["secs"] > 2400 and x["totM"] < 40_000_000 for x in runs[7]), "d=7 runs stopped by the time limit"
put("m.d7.timelimited", len(runs[7]), SRC, "d=7 runs that ended at the 2400 s limit with se > 2% (all of them)")
put("m.lowd.runs", len(low), SRC, "d=3,5 runs, all stopped by the 2% rule")
put("m.lam.max", f"{max(x['lam'] for v in runs.values() for x in v):.1f}", SRC, "largest mean number of faults sum q_i over the 12 cells (d=7, p=1e-3)")
H = [json.load(open(f)) for f in sorted(glob.glob(f"{R}/holdout_truth_d9_p1e-4_s*.json"))]; assert len(H) == 6
put("m.N.d9", H[0]["N"], "results/holdout_truth_d9_p1e-4_s*.json", "fault mechanisms, d=9")
assert all(x["strata"][str(w)][3] == 0 and x["strata"][str(w)][2] == 100000 for x in H for w in range(1, 5))
put("m.h.secs.min", f"{min(x['secs'] for x in H):.0f}", "results/holdout_truth_d9_p1e-4_s*.json", "holdout run time, s")
put("m.h.totM.min", f"{min(x['totM'] for x in H)/1e6:.1f}", "results/holdout_truth_d9_p1e-4_s*.json", "holdout samples per run, millions")
put("m.h.totM.max", f"{max(x['totM'] for x in H)/1e6:.1f}", "results/holdout_truth_d9_p1e-4_s*.json", "")
# certified-interval check on the surrogate (registered: PLAN.md 'Certified interval amendment'); truth = exact polynomial
import sys; sys.path.insert(0, "../code")
CS = "results/cert_validation/*/cert_exact_*.json vs exact/d{3,5}_manifest.json"
bs = []; wid = []; tot_in = tot_n = 0
for d in (3, 5):
    c = json.load(open(f"../exact/d{d}_manifest.json"))["coefficients_p_L=sum c_k p^k"]
    for p, t in ((1e-3, "3"), (1e-4, "4")):
        ex = sum(v["coef"] * p ** int(k) for k, v in c.items())
        rs = [json.load(open(f)) for f in sorted(glob.glob(f"{R}/cert_validation/d{d}_p1e-{t}/cert_exact_d{d}_p{p:g}_r*.json"))]
        assert len(rs) == 24 and all(r["kind"] == "exact" and r["delta"] == 0.05 for r in rs)
        ins = sum(r["L"] <= ex <= r["U"] for r in rs); rw = float(np.median([(r["U"] - r["L"]) / ex for r in rs]))
        put(f"c.d{d}.t{t}.inside", ins, CS, "replicates with exact value inside [L,U]"); put(f"c.d{d}.t{t}.n", len(rs), CS, "")
        put(f"c.d{d}.t{t}.width", f"{100*rw:.1f}", CS, "median (U-L)/exact, percent"); wid.append(rw); tot_in += ins; tot_n += len(rs)
        bs += [s["b"] for r in rs for s in r["per"].values()]
put("c.inside.all", tot_in, CS, ""); put("c.n.all", tot_n, CS, "")
put("c.width.min", f"{100*min(wid):.1f}", CS, ""); put("c.width.max", f"{100*max(wid):.1f}", CS, "")
A = json.load(open(f"{R}/cert_d9_alloc.json")); bs9 = list(A["b"].values())
put("c.b.min", f"{min(bs + bs9):.2f}", CS + " + results/cert_d9_alloc.json", "smallest rigorous weight bound b_w over all strata used"); put("c.b.max", f"{max(bs + bs9):.2f}", CS + " + results/cert_d9_alloc.json", "largest")
put("c.d9.wmax", A["wmax"], "results/cert_d9_alloc.json", "largest stratum bounded at d=9"); put("c.d9.ntot", f"10^{{{int(round(math.log10(A['Ntot'])))}}}", "results/cert_d9_alloc.json", "registered sample budget")
# (v5) the d=9 certified result now exists and is reported from make_numbers_v5.py; the v3 guard on its absence is removed.
# (v4) the below-crossing results now exist and are reported from make_numbers_v4.py; the v3 guard on their absence is removed.
with open("numbers_v3.tex", "w") as f:
    f.write("% GENERATED by make_numbers_v3.py. Do not edit. Input after numbers_amend.tex.\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers_v3.json", "w"), indent=1)
for k, (v, s, h) in V.items(): print(k, v)
