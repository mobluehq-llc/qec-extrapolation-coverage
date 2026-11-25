"""Every number in the paper, from the committed result files. Writes numbers.tex (macros), tab_*.tex, CLAIMS.md, numbers.json.
Run from this folder:  python3 make_numbers.py
Nothing is simulated here. The frozen analysis code (../code/coverage.py) is imported unchanged; the coverage CSVs written on a compute host are
re-derived on this machine and asserted equal (hits per cell) before any number is emitted."""
import sys, json, csv, glob, re, math, os
import numpy as np
sys.path.insert(0, "../code")
from coverage import fit_all, load_reps, wls, cp, PS, X0, TARGETS, coverage  # frozen code

R = "../results"; EX = "../exact"
V = {}      # key -> (tex string, source, how)
def put(key, val, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", key), key
    assert key not in V, key
    V[key] = (str(val), src, how)
def sci(x, nd=2):
    if x == 0: return "0"
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
def pct(x, nd=1): return f"{100*x:.{nd}f}"
def fx(x):  # ratio formatting: 2 significant figures below 10, integer above
    return f"{x:.2f}" if x < 1 else (f"{x:.1f}" if x < 10 else f"{x:.0f}")
RC = ["R1", "R2", "R3", "R4", "R5"]; DS = [3, 5, 7]; TG = [1e-3, 1e-4, 1e-5, 1e-6]
tn = {1e-3: "3", 1e-4: "4", 1e-5: "5", 1e-6: "6"}     # target label: decade exponent

# ---------------------------------------------------------------- truth (12 cells), as analyze_main.load_truth
T = {}
for f in glob.glob(f"{R}/truth_d*_p*.json"):
    x = json.load(open(f)); T[(x["d"], x["p"])] = x
for f in glob.glob(f"{R}/truthextra_d*_p*.json"):
    y = json.load(open(f)); x = T[(y["d"], y["p"])]; w1, w2 = 1 / x["se"] ** 2, 1 / y["se"] ** 2
    x["t1"], x["s1"], x["t2"], x["s2"] = x["truth"], x["se"], y["truth"], y["se"]
    x["truth"] = (w1 * x["t1"] + w2 * y["truth"]) / (w1 + w2); x["se"] = (w1 + w2) ** -0.5
truth = {(d, p): T[(d, p)]["truth"] for d in DS for p in TG}
relse = {(d, p): T[(d, p)]["se"] / T[(d, p)]["truth"] for d in DS for p in TG}
zmax = 0
for d in DS:
    for p in TG:
        x = T[(d, p)]; z = (x["t1"] - x["t2"]) / math.hypot(x["s1"], x["s2"]); zmax = max(zmax, abs(z))
        put(f"truth.d{d}.t{tn[p]}", sci(x["truth"], 3), f"results/truth_d{d}_p*.json + truthextra", "inverse-variance combination of run 1 and run 2")
        put(f"truthse.d{d}.t{tn[p]}", pct(relse[(d, p)]), "same", "combined se / combined estimate, percent")
put("truth.zmax", f"{zmax:.2f}", "results/truth_*, truthextra_*", "max |z| between run 1 and run 2 over the 12 cells")
put("truthse.lowd.min", pct(min(relse[(d, p)] for d in (3, 5) for p in TG)), "same", "min rel se, d=3,5")
put("truthse.lowd.max", pct(max(relse[(d, p)] for d in (3, 5) for p in TG)), "same", "max rel se, d=3,5")
put("truthse.d7.far.min", pct(min(relse[(7, p)] for p in TG[1:]), 0), "same", "min rel se, d=7, targets 1e-4..1e-6")
put("truthse.d7.far.max", pct(max(relse[(7, p)] for p in TG[1:]), 0), "same", "max rel se, d=7, targets 1e-4..1e-6")

# ---------------------------------------------------------------- truth validation against direct Monte Carlo
rows = list(csv.DictReader(open(f"{R}/L_truth_validation.csv")))
el = [r for r in rows if r["eligible"] == "1"]
put("val.rows", len(rows), "results/L_truth_validation.csv", "data rows")
put("val.eligible", len(el), "same", "rows with eligible=1")
put("val.ok", sum(r["ok"] == "1" for r in el), "same", "eligible rows with |z|<=2")
put("val.zmax", f"{max(abs(float(r['z'])) for r in el):.2f}", "same", "max |z|")
put("val.minerr", min(int(r["direct_errors"]) for r in el), "same", "fewest direct-MC logical errors in a cell")
put("val.pmin", sci(min(float(r["truth"]) for r in el), 2), "same", "smallest validated p_L")

# ---------------------------------------------------------------- coverage (frozen CSVs; re-derived here and compared)
reps = load_reps(f"{R}/fit_d357.json")
put("reps", len(reps), "results/fit_d357.json", "replicates")
put("windowpoints", sum(1 for _ in json.load(open(f"{R}/fit_d357.json"))), "same", "(rep,d,window p) points")
COV = {}
for scaled, tag, fn in ((False, "pri", "L_coverage.csv"), (True, "sca", "L_coverage_scaled.csv")):
    cs = {(r["recipe"], int(r["d"]), float(r["target"])): r for r in csv.DictReader(open(f"{R}/{fn}"))}
    mine = {(r["recipe"], r["d"], r["target"]): r for r in coverage(reps, truth, scaled=scaled)}
    assert set(cs) == set(mine) and all(int(cs[k]["hits"]) == mine[k]["hits"] for k in cs), f"re-derived coverage differs from {fn}"
    COV[tag] = cs
    for (rc, d, p), r in cs.items():
        put(f"cov.{tag}.{rc}.d{d}.t{tn[p]}", pct(float(r["cov"])), f"results/{fn}", "hits/300, percent")
    c = [float(r["cov"]) for r in cs.values()]
    put(f"cells.{tag}", len(c), f"results/{fn}", "cells")
    put(f"cells.{tag}.zero", sum(x == 0 for x in c), f"results/{fn}", "cells with 0 hits")
    put(f"cells.{tag}.below90", sum(x < 0.9 for x in c), f"results/{fn}", "cells with coverage < 0.90")
    put(f"cells.{tag}.firm", sum(float(r["cp_hi"]) < 0.9 for r in cs.values()), f"results/{fn}", "cells with CP upper < 0.90")
    put(f"cells.{tag}.calibrated", sum(0.93 <= x <= 0.97 for x in c), f"results/{fn}", "cells in 93-97%")
    put(f"cells.{tag}.above97", sum(x > 0.97 for x in c), f"results/{fn}", "cells above 97%")
    put(f"cov.{tag}.max", pct(max(c)), f"results/{fn}", "largest coverage")
    for rc in RC:
        put(f"cov.{tag}.median.{rc}", pct(float(np.median([float(r["cov"]) for k, r in cs.items() if k[0] == rc]))), f"results/{fn}", "median over the recipe's cells")
put("cp.zero.upper", pct(cp(0, 300)[1]), "code/coverage.py cp(0,300)", "Clopper-Pearson upper limit for 0/300")
far = TG[1:]
for rc in ("R1", "R5"):
    c = [r for k, r in COV["pri"].items() if k[0] == rc and k[2] in far]
    put(f"h1.{rc}.cells", len(c), "results/L_coverage.csv", "cells with target <= 1e-4")
    put(f"h1.{rc}.mis", sum(float(r["cov"]) < 0.9 for r in c), "same", "of which coverage < 0.90")
    put(f"h1.{rc}.firm", sum(float(r["cp_hi"]) < 0.9 for r in c), "same", "of which CP upper < 0.90")
vt = open(f"{R}/L_verdicts.txt").read()
assert "primary(absolute sigma) H1 (R1 and R5 miscalibrated in >=2/3 of cells at >=2 decades): SUPPORTED" in vt
assert "sensitivity(chi2/dof-scaled) H1 (R1 and R5 miscalibrated in >=2/3 of cells at >=2 decades): SUPPORTED" in vt

# ---------------------------------------------------------------- direction and interval width (same frozen fits)
LR = {}; HW = {}
for data in reps:
    for p in TG:
        o = fit_all({d: data[d] for d in DS}, p)
        for (rc, d), (mu, var) in o.items():
            LR.setdefault((rc, d, p), []).append(mu - math.log(truth[(d, p)])); HW.setdefault((rc, d, p), []).append(1.96 * math.sqrt(var))
RAT = {k: float(np.exp(np.median(v))) for k, v in LR.items()}
for (rc, d, p), v in RAT.items():
    put(f"ratio.{rc}.d{d}.t{tn[p]}", fx(v), "results/fit_d357.json, truth files; code/coverage.py fit_all", "median over 300 replicates of prediction/truth")
for rc in ("R2", "R3"):
    v = [RAT[(rc, d, p)] for d in DS for p in TG]
    put(f"ratio.{rc}.min", f"{min(v):.2f}", "same", "min over 12 cells"); put(f"ratio.{rc}.max", f"{max(v):.2f}", "same", "max over 12 cells")
v = [RAT[(rc, d, p)] for rc in ("R2", "R3") for d in DS for p in TG]
put("ratio.fixed.min", f"{min(v):.2f}", "same", "min over R2,R3 cells"); put("ratio.fixed.max", f"{max(v):.2f}", "same", "max over R2,R3 cells")
put("ratio.fixed.lowpct.max", f"{100*(1-min(v)):.0f}", "same", "largest shortfall, percent"); put("ratio.fixed.lowpct.min", f"{100*(1-max(v)):.0f}", "same", "smallest shortfall, percent")
put("ratio.R4.min", fx(min(RAT[("R4", d, 1e-3)] for d in DS)), "same", "R4 at 1e-3, min over d"); put("ratio.R4.max", fx(max(RAT[("R4", d, 1e-6)] for d in DS)), "same", "R4 at 1e-6, max over d")
put("ratio.max", fx(max(RAT.values())), "same", "largest median ratio in the table")
for rc in RC:
    hw = [float(np.median(HW[k])) for k in HW if k[0] == rc]
    put(f"hw.{rc}.min", pct(math.exp(min(hw)) - 1), "same", "smallest cell-median 95% half-width, percent of prediction")
    put(f"hw.{rc}.max", pct(math.exp(max(hw)) - 1, 0 if math.exp(max(hw)) - 1 > 0.2 else 1), "same", "largest cell-median 95% half-width, percent of prediction")
# sigma distance of the fixed-exponent forms
zz = [abs(np.median(LR[(rc, d, p)])) / (np.median(HW[(rc, d, p)]) / 1.96) for rc in ("R2", "R3") for d in (3, 5) for p in TG]
put("fixed.sigmas.min", f"{min(zz):.0f}", "same", "R2,R3 at d=3,5: |median log error| / median interval sigma, min"); put("fixed.sigmas.max", f"{max(zz):.0f}", "same", "max")

# descriptive: window values, crossing, fitted R1 slope, chi2/dof
for d in DS:
    P = np.array([r[d][0] / r[d][1] for r in reps]).mean(0)
    for ip in range(5): put(f"win.d{d}.p{ip}", sci(P[ip], 3), "results/fit_d357.json", "mean over replicates of k/n")
    sl = []; c2 = []
    for r in reps:
        k, n = r[d]; Pq = k / n; b, C, chi2, dof = wls(np.c_[np.ones(5), X0], np.log(Pq), k / (1 - Pq)); sl.append(b[1]); c2.append(chi2 / dof)
    put(f"slope.d{d}", f"{np.median(sl):.2f}", "same; code/coverage.py wls", "median fitted R1 exponent"); put(f"chidof.d{d}", f"{np.median(c2):.1f}", "same", "median chi2/dof of the R1 fit")
    put(f"nominal.d{d}", f"{(d+1)/2:.0f}", "definition", "(d+1)/2")
Pm = {d: np.array([r[d][0] / r[d][1] for r in reps]).mean(0) for d in DS}
order = ["".join(str(d) for d in sorted(DS, key=lambda d: Pm[d][ip])) for ip in range(5)]   # ascending p_L
assert order[0] == "753" and order[1] == "753" and order[3] == "357" and order[4] == "357", order
put("win.order.p2", order[2], "results/fit_d357.json", "distances in ascending order of mean p_L at p=8e-3")

# EXPLORATORY (post hoc, not registered): propagate the truth se into the hit test
hits = {}
for data in reps:
    for p in TG:
        o = fit_all({d: data[d] for d in DS}, p)
        for (rc, d), (mu, var) in o.items():
            hits[(rc, d, p)] = hits.get((rc, d, p), 0) + int(abs(math.log(truth[(d, p)]) - mu) <= 1.96 * math.sqrt(var + relse[(d, p)] ** 2))
c = {k: v / 300 for k, v in hits.items()}
put("x.prop.below90", sum(x < 0.9 for x in c.values()), "EXPLORATORY; results/fit_d357.json, truth files", "cells < 0.90 when truth se is added in quadrature")
put("x.prop.zero", sum(x == 0 for x in c.values()), "same", "cells at 0")
put("x.prop.changed", sum(abs(c[k] - float(COV["pri"][k]["cov"])) > 0.05 for k in c), "same", "cells whose coverage moves by more than 5 points")
ch = sorted([k for k in c if c[k] >= 0.9]); put("x.prop.above90.list", ", ".join(f"{k[0]} $d={k[1]}$ $10^{{-{tn[k[2]]}}}$" for k in ch) or "none", "same", "cells reaching 0.90")
for k in c:
    put(f"x.prop.{k[0]}.d{k[1]}.t{tn[k[2]]}", pct(c[k]), "same", "coverage with truth se propagated, percent")
put("x.prop.lowd.changed", sum(abs(c[k] - float(COV["pri"][k]["cov"])) > 0.05 for k in c if k[1] != 7), "same", "d=3,5 cells moving by more than 5 points")

# EXPLORATORY (post hoc): refit R1 and R2 on the three lowest window points only (4e-3, 6e-3, 8e-3)
def fit_low(data, tp, m=3):
    out = {}; xs = math.log(tp)
    for d in data:
        k, n = data[d]; k = k[:m]; n = n[:m]; P = k / n; y = np.log(P); w = k / (1 - P); x0 = X0[:m]
        b, C, c2, dof = wls(np.c_[np.ones(m), x0], y, w); g = np.array([1, xs]); out[("R1", d)] = (g @ b, g @ C @ g, b[1])
        h = (d + 1) / 2; b, C, c2, dof = wls(np.ones((m, 1)), y - h * x0, w); out[("R2", d)] = (b[0] + h * xs, C[0, 0], h)
    return out
lh = {}; lr = {}; ls = {}
for data in reps:
    for p in TG:
        for (rc, d), (mu, var, s) in fit_low({d: data[d] for d in DS}, p).items():
            lh[(rc, d, p)] = lh.get((rc, d, p), 0) + int(abs(math.log(truth[(d, p)]) - mu) <= 1.96 * math.sqrt(var)); lr.setdefault((rc, d, p), []).append(mu - math.log(truth[(d, p)])); ls.setdefault((rc, d), []).append(s)
for k in lh:
    put(f"x.low.cov.{k[0]}.d{k[1]}.t{tn[k[2]]}", pct(lh[k] / 300), "EXPLORATORY; results/fit_d357.json (three lowest window points), truth files", "coverage, percent")
    put(f"x.low.ratio.{k[0]}.d{k[1]}.t{tn[k[2]]}", fx(float(np.exp(np.median(lr[k])))), "same", "median prediction/truth")
put("x.low.below90", sum(v / 300 < 0.9 for v in lh.values()), "same", "cells below 0.90"); put("x.low.cells", len(lh), "same", "cells")
put("x.low.cov.max", pct(max(lh.values()) / 300), "same", "largest coverage")
for d in DS: put(f"x.low.slope.d{d}", f"{np.median(ls[('R1', d)]):.2f}", "same", "median R1 exponent on three lowest points")
put("x.low.R1.max", fx(max(float(np.exp(np.median(lr[k]))) for k in lr if k[0] == "R1")), "same", "largest R1 median ratio")
put("x.low.R2.min", fx(min(float(np.exp(np.median(lr[k]))) for k in lr if k[0] == "R2")), "same", "smallest R2 median ratio")
put("x.low.R2.max", fx(max(float(np.exp(np.median(lr[k]))) for k in lr if k[0] == "R2")), "same", "largest R2 median ratio")

# ---------------------------------------------------------------- controls (text written by the frozen controls.py)
ta = open(f"{R}/L_control_a.txt").read()
for rc in ("R1", "R2", "R3", "R4"):
    m = re.search(rf"CONTROL-A {rc} dataonly: reject rate ([\d.]+); cells with reject rate>=0.95: (\d+)/(\d+)", ta)
    put(f"ca.{rc}.do", m.group(1), "results/L_control_a.txt", "reject rate, data-only dataset"); put(f"ca.{rc}.do.cells", m.group(2), "same", "cells >= 0.95"); assert m.group(3) == "15"
    m = re.search(rf"CONTROL-A {rc} circuit_fresh: reject rate ([\d.]+); cells with reject rate>=0.95: (\d+)/(\d+)", ta)
    put(f"ca.{rc}.fresh", pct(float(m.group(1))), "same", "reject rate, fresh circuit data, percent"); put(f"ca.{rc}.fresh.cells", m.group(2), "same", "cells >= 0.95")
tb = open(f"{R}/L_control_b.txt").read()
for rc in RC:
    m = re.search(rf"CONTROL-B {rc}: miscalibrated (\d+)/(\d+) cells; median coverage ([\d.]+)", tb)
    put(f"cb.{rc}.mis", m.group(1), "results/L_control_b.txt", "miscalibrated cells"); put(f"cb.{rc}.cells", m.group(2), "same", "cells"); put(f"cb.{rc}.median", pct(float(m.group(3))), "same", "median coverage, percent")

# ---------------------------------------------------------------- exact enumeration (linearised-DEM surrogate)
M = {3: json.load(open(f"{EX}/d3_manifest.json")), 5: json.load(open(f"{EX}/d5_manifest.json")), 7: json.load(open(f"{EX}/d7_manifest.json"))}
S7 = json.load(open(glob.glob(f"{EX}/d7_w4_snapshot_*.json")[0]))
def big(n): return f"{n:,}".replace(",", "{,}")
for d in DS:
    m = M[d]; put(f"ex.d{d}.nmech", big(m["n_mechanisms"]), f"exact/d{d}_manifest.json", "mechanisms"); put(f"ex.d{d}.ndet", m["n_detectors"], f"exact/d{d}_manifest.json", "detectors")
    for w, s in m["shells"].items():
        assert s["status"] == "COMPLETE" and s["sets_done"] == s["sets_total"]
        put(f"ex.d{d}.w{w}.sets", big(s["sets_total"]) if s["sets_total"] < 1e9 else sci(s["sets_total"], 3), f"exact/d{d}_manifest.json", "fault sets in the shell (all enumerated)")
        put(f"ex.d{d}.w{w}.fail", big(s["n_fail"]) if s["n_fail"] < 1e9 else sci(s["n_fail"], 3), f"exact/d{d}_manifest.json", "failing sets")
        put(f"ex.d{d}.w{w}.frac", "0" if s["n_fail"] == 0 else sci(s["n_fail"] / s["sets_total"], 2), f"exact/d{d}_manifest.json", "failing fraction")
    for k, cdict in m["coefficients_p_L=sum c_k p^k"].items():
        assert cdict["status"] == "EXACT"
        cv = cdict["coef"]; put(f"ex.d{d}.c{k}", "0" if cv == 0 else (f"{cv:.2f}" if abs(cv) < 1e4 else f"{cv:.0f}").replace("-", "-"), f"exact/d{d}_manifest.json", "exact coefficient")
c3, c4 = (M[5]["coefficients_p_L=sum c_k p^k"][k]["coef"] for k in ("3", "4"))
put("ex.d5.c4overc3", f"{100*abs(c4)*1e-3/c3:.2f}", "exact/d5_manifest.json", "|c4| p / c3 at p=1e-3, percent")
d3c = {int(k): v["coef"] for k, v in M[3]["coefficients_p_L=sum c_k p^k"].items()}
put("ex.d3.c3overc2", f"{100*abs(d3c[3])*1e-3/d3c[2]:.2f}", "exact/d3_manifest.json", "|c3| p / c2 at p=1e-3, percent")
put("ex.d7.w4.units", S7["units_done"], "exact/d7_w4_snapshot_*.json", "sampled units decoded"); put("ex.d7.w4.unitstotal", big(S7["units_total"]), "same", "units in the shell")
put("ex.d7.w4.setsdone", sci(S7["sets_done"], 2), "same", "sets decoded"); put("ex.d7.w4.setstotal", sci(S7["sets_total"], 2), "same", "sets in the shell")
put("ex.d7.w4.fracdone", f"{100*S7['frac_sets_done']:.2f}", "same", "percent of the shell decoded")
put("ex.d7.w4.fail", big(S7["observed"]["n_fail"]), "same", "failing sets observed")
put("ex.d7.w4.frac", sci(S7["observed"]["n_fail"] / S7["sets_done"], 2), "same", "observed failing fraction (sets done)")
e4 = S7["ESTIMATE_full_shell_NOT_exact"]["T0"]
put("ex.d7.c4.est", sci(e4["estimate"], 2), "same", "design-based estimate of c4"); put("ex.d7.c4.se", sci(e4["se"], 2), "same", "its se"); put("ex.d7.c4.relse", f"{100*e4['rel_se']:.0f}", "same", "rel se, percent")
put("ex.d7.c4.lower", f"{S7['observed']['T'][0]:.0f}", "same", "observed lower bound on c4")
put("ex.d7.c4.circuit", sci(truth[(7, 1e-6)] / 1e-24, 2), "results/truth_d7_p1e-6.json + truthextra", "real-circuit truth at p=1e-6 divided by p^4")
def poly(d, p): return sum(v["coef"] * p ** int(k) for k, v in M[d]["coefficients_p_L=sum c_k p^k"].items())
rr = []
for d in (3, 5):
    for p in TG:
        r = truth[(d, p)] / poly(d, p); rr.append((d, p, r))
        put(f"ex.gap.d{d}.t{tn[p]}", f"{100*(r-1):+.1f}", f"truth files; exact/d{d}_manifest.json", "real-circuit truth / surrogate polynomial - 1, percent")
        put(f"ex.poly.d{d}.t{tn[p]}", sci(poly(d, p), 3), f"exact/d{d}_manifest.json", "polynomial through order 4")
put("ex.gap.d5.mean", f"{100*(np.mean([r for d, p, r in rr if d == 5])-1):.0f}", "same", "mean over 4 targets, d=5, percent")
put("ex.gap.d3.mean", f"{100*(np.mean([r for d, p, r in rr if d == 3])-1):.0f}", "same", "mean over 4 targets, d=3, percent")

# ---------------------------------------------------------------- IS validation against the exact polynomial
IS = {}
for d, p, src in ((3, "1e-3", "d3_p1e-3.log"), (3, "1e-4", "d3_p1e-4.log")):
    IS[(d, p)] = json.loads(open(f"{R}/is_validation/{src}").read().strip().splitlines()[-1]); IS[(d, p)]["src"] = f"results/is_validation/{src}"
for line in open(f"{R}/is_validation/d5_agg.json"):
    if line.strip():
        x = json.loads(line); IS[(5, f"{x['p']:.0e}".replace("e-0", "e-"))] = x; x["src"] = "results/is_validation/d5_agg.json"
assert len(IS) == 4
for (d, p), x in IS.items():
    t = p[-1]; assert x["verdict"] == "PASS"
    put(f"is.d{d}.t{t}.exact", sci(x["exact"], 5), x["src"], "exact polynomial value"); put(f"is.d{d}.t{t}.inside", x["inside"], x["src"], "seeds with the exact value inside the 95% interval")
    put(f"is.d{d}.t{t}.n", x["n_runs"], x["src"], "seeds"); put(f"is.d{d}.t{t}.width", pct(x["median_rel_width"]), x["src"], "median relative width, percent")
    assert abs(x["exact"] / poly(d, float(p)) - 1) < 1e-9
put("is.width.min", pct(min(x["median_rel_width"] for x in IS.values())), "results/is_validation", "min of the four medians"); put("is.width.max", pct(max(x["median_rel_width"] for x in IS.values())), "same", "max")

# ---------------------------------------------------------------- d=9 holdout
H = [json.load(open(f)) for f in sorted(glob.glob(f"{R}/holdout_truth_d9_p1e-4_s*.json"))]
assert len(H) == 6
w = np.array([1 / r["se"] ** 2 for r in H]); t = np.array([r["truth"] for r in H]); ht = float((w * t).sum() / w.sum()); hs = float(w.sum() ** -0.5)
put("h.runs", len(H), "results/holdout_truth_d9_p1e-4_s900*.json", "truth runs"); put("h.truth", sci(ht, 2), "same; code/holdout.py", "inverse-variance combination"); put("h.se", sci(hs, 2), "same", "combined se")
put("h.relse", f"{hs/ht:.0f}", "same", "combined se / estimate (a multiple, not percent)")
ts = np.sort(t); put("h.five.min", sci(ts[0], 2), "same", "smallest run estimate"); put("h.five.max", sci(ts[4], 2), "same", "fifth smallest run estimate"); put("h.outlier", sci(ts[5], 2), "same", "largest run estimate")
rs = np.array([r["se"] / r["truth"] for r in H]); put("h.run.relse.min", f"{rs.min():.0f}", "same", "min per-run se/estimate"); put("h.run.relse.max", f"{rs.max():.0f}", "same", "max per-run se/estimate")
put("h.samples.min", f"{min(r['totM'] for r in H)/1e6:.1f}", "same", "samples per run, millions, min"); put("h.samples.max", f"{max(r['totM'] for r in H)/1e6:.1f}", "same", "max")
put("h.upper", sci(ht + 1.96 * hs, 2), "same", "estimate + 1.96 se")
hc = {r["recipe"]: r for r in csv.DictReader(open(f"{R}/L_holdout_d9.csv"))}
for rc in RC: put(f"h.hits.{rc}", hc[rc]["hits"], "results/L_holdout_d9.csv", "hits of 300"); assert hc[rc]["reps"] == "300"
r9 = load_reps(f"{R}/fit_d9.json"); allr = [{**a, **b} for a, b in zip(reps, r9)]
hm = {rc: [] for rc in RC}
for data in allr:
    o = fit_all({d: data[d] for d in (3, 5, 9)}, 1e-4)
    for rc in RC: hm[rc].append(o[(rc, 9)][0])
for rc in RC:
    med = float(np.exp(np.median(hm[rc]))); put(f"h.pred.{rc}", sci(med, 2), "results/fit_d357.json, fit_d9.json; code/coverage.py", "median prediction of p_L(d=9, 1e-4)")
    put(f"h.ratio.{rc}", fx(med / ht), "same + holdout truth", "median prediction / point estimate of truth")
    put(f"h.above.{rc}", sum(math.exp(x) > ht + 1.96 * hs for x in hm[rc]), "same", "replicates (of 300) whose prediction exceeds estimate + 1.96 se")
    put(f"h.aboveout.{rc}", sum(math.exp(x) > ts[5] for x in hm[rc]), "same", "replicates whose prediction exceeds the largest single-run estimate")

# ---------------------------------------------------------------- outputs
with open("numbers.tex", "w") as f:
    f.write("% GENERATED by make_numbers.py from the committed result files. Do not edit.\n")
    f.write("\\makeatletter\\newcommand{\\V}[1]{\\@ifundefined{LV@#1}{\\errmessage{undefined number #1}}{\\@nameuse{LV@#1}}}\\makeatother\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers.json", "w"), indent=1)

def covtab(tag, fn):
    L = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
    for rc in RC:
        for d in DS:
            if rc == "R5" and d != 7: continue
            L.append(f"{rc} & {d} & " + " & ".join(f"\\V{{cov.{tag}.{rc}.d{d}.t{tn[p]}}}" for p in TG) + "\\\\")
    L += ["\\hline\\hline", "\\end{tabular}"]; open(fn, "w").write("\n".join(L) + "\n")
covtab("pri", "tab_cov_primary.tex"); covtab("sca", "tab_cov_scaled.tex")
L = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
for rc in RC:
    for d in DS:
        if rc == "R5" and d != 7: continue
        L.append(f"{rc} & {d} & " + " & ".join(f"\\V{{ratio.{rc}.d{d}.t{tn[p]}}}" for p in TG) + "\\\\")
L += ["\\hline\\hline", "\\end{tabular}"]; open("tab_ratio.tex", "w").write("\n".join(L) + "\n")
L = ["\\begin{tabular}{lrrrr}\\hline\\hline", "$d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
for d in DS: L.append(f"{d} & " + " & ".join(f"$\\V{{truth.d{d}.t{tn[p]}}}$ (\\V{{truthse.d{d}.t{tn[p]}}})" for p in TG) + "\\\\")
L += ["\\hline\\hline", "\\end{tabular}"]; open("tab_truth.tex", "w").write("\n".join(L) + "\n")
L = ["\\begin{tabular}{llrrr}\\hline\\hline", "$d$ & weight & fault sets & failing sets & failing fraction\\\\\\hline"]
for d in DS:
    for w in M[d]["shells"]: L.append(f"{d} & {w} & $\\V{{ex.d{d}.w{w}.sets}}$ & $\\V{{ex.d{d}.w{w}.fail}}$ & $\\V{{ex.d{d}.w{w}.frac}}$\\\\")
L.append("7 & 4 (sampled) & $\\V{ex.d7.w4.setsdone}$ of $\\V{ex.d7.w4.setstotal}$ & $\\V{ex.d7.w4.fail}$ & $\\V{ex.d7.w4.frac}$\\\\")
L += ["\\hline\\hline", "\\end{tabular}"]; open("tab_shells.tex", "w").write("\n".join(L) + "\n")
for name, pre in (("tab_x_prop.tex", "x.prop."),):
    L = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
    for rc in RC:
        for d in DS:
            if rc == "R5" and d != 7: continue
            L.append(f"{rc} & {d} & " + " & ".join(f"\\V{{{pre}{rc}.d{d}.t{tn[p]}}}" for p in TG) + "\\\\")
    L += ["\\hline\\hline", "\\end{tabular}"]; open(name, "w").write("\n".join(L) + "\n")
L = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
for rc in ("R1", "R2"):
    for d in DS: L.append(f"{rc} & {d} & " + " & ".join(f"\\V{{x.low.cov.{rc}.d{d}.t{tn[p]}}} (\\V{{x.low.ratio.{rc}.d{d}.t{tn[p]}}})" for p in TG) + "\\\\")
L += ["\\hline\\hline", "\\end{tabular}"]; open("tab_x_low.tex", "w").write("\n".join(L) + "\n")
L = ["\\begin{tabular}{lrrrrr}\\hline\\hline", "$d$ & $4\\times10^{-3}$ & $6\\times10^{-3}$ & $8\\times10^{-3}$ & $10^{-2}$ & $1.5\\times10^{-2}$\\\\\\hline"]
for d in DS: L.append(f"{d} & " + " & ".join(f"$\\V{{win.d{d}.p{ip}}}$" for ip in range(5)) + "\\\\")
L += ["\\hline\\hline", "\\end{tabular}"]; open("tab_window.tex", "w").write("\n".join(L) + "\n")

with open("CLAIMS.md", "w") as f:
    f.write("# Paper L: number-to-source list (GENERATED by paper/make_numbers.py; do not edit)\n\n")
    f.write("Every number in `main.tex` and `supplement.tex` is typeset through `\\V{key}`; an undefined key stops the build. Paths are relative to `ideas/L_qec_extrapolation_calibration/`. Rows tagged EXPLORATORY are post hoc and were not registered. The frozen claims ledger is `../CLAIMS.md` and is not changed by this file.\n\n")
    f.write("**Discrepancy found while building this list:** `RESULTS.md` and `../CLAIMS.md` state the truth validation as 17/17 cells; `results/L_truth_validation.csv` has 17 lines including its header, i.e. 16 data rows, all eligible and all with |z| <= 2. The paper reports 16/16.\n\n")
    f.write("| key | value | source | how |\n|---|---|---|---|\n")
    for k, (v, s, h) in V.items(): f.write(f"| `{k}` | {v.replace('|','/')} | {s} | {h} |\n")
print(len(V), "numbers written")
for k in ["cells.pri.zero", "cells.pri.below90", "cov.pri.max", "ratio.max", "ratio.fixed.min", "ratio.fixed.max", "hw.R1.min", "hw.R1.max", "hw.R2.min", "hw.R2.max", "hw.R3.min", "hw.R3.max", "hw.R4.max", "hw.R5.max", "fixed.sigmas.min", "fixed.sigmas.max",
          "x.prop.below90", "x.prop.zero", "x.prop.changed", "x.prop.above90.list", "x.prop.lowd.changed", "x.low.below90", "x.low.cov.max", "x.low.R1.max", "x.low.R2.min", "x.low.R2.max", "x.low.slope.d3", "x.low.slope.d5", "x.low.slope.d7",
          "ex.d5.c4overc3", "ex.d3.c3overc2", "ex.d7.c4.est", "ex.d7.c4.circuit", "ex.gap.d5.mean", "ex.gap.d3.mean", "h.truth", "h.se", "h.relse", "h.upper", "h.five.min", "h.five.max", "h.outlier", "h.run.relse.min", "h.run.relse.max",
          "truth.zmax", "val.eligible", "val.ok", "val.zmax", "win.order.p2", "truthse.lowd.min", "truthse.lowd.max", "truthse.d7.far.min", "truthse.d7.far.max", "truthse.d7.t3"] + [f"h.ratio.{r}" for r in RC] + [f"h.pred.{r}" for r in RC] + [f"h.above.{r}" for r in RC] + [f"h.aboveout.{r}" for r in RC] + [f"ex.gap.d{d}.t{t}" for d in (3, 5) for t in "3456"] + [f"slope.d{d}" for d in DS] + [f"chidof.d{d}" for d in DS]:
    print(k, "=", V[k][0])
