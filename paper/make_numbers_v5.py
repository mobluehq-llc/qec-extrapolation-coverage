"""Numbers for draft v5: (1) the certified interval at d=9, p=1e-4 (PLAN.md 'CERTIFIED INTERVAL AMENDMENT'; commit d829cc8),
(2) R4' and R2-lin on the window below the crossing (PLAN.md 'AMENDMENT R4P-BW'; commit f2619b8),
(3) rolling-origin discrepancy inflation for R2/R3 (PLAN.md 'AMENDMENT 2'; commit a438b4c).
Run from this folder after make_numbers.py, make_numbers_amend.py, make_numbers_v3.py and make_numbers_v4.py:
  python3 make_numbers_v5.py
Reads only stored result files. Nothing is simulated and no stored result is changed. The counts quoted in RESULTS.md,
RESULTS_R4P_BW.md and RESULTS_ROLLING_INFLATION.md are asserted here before any number is emitted.
Writes numbers_v5.tex (keys c9.*, r4b.*, inf.*), numbers_v5.json, tab_r4p_bw.tex, tab_cert_d9.tex, tab_cert_recipes.tex, tab_inflate.tex."""
import sys, json, csv, glob, math, re, os
import numpy as np
sys.path.insert(0, "../code")
os.chdir("../code")
from coverage import fit_all, load_reps
os.chdir("../paper")
R = "../results"; V = {}
def put(k, v, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), src, how)
def pct(x, nd=1): return f"{100*x:.{nd}f}"
def sci(x, nd=2):
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
def sci_out(x, up, nd=3):
    """round OUTWARD to nd significant figures: up=True ceil (upper endpoint), False floor (lower endpoint)"""
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** (e - nd + 1); r = math.ceil(m - 1e-9) if up else math.floor(m + 1e-9)
    if r >= 10 ** nd: r = r // 10; e += 1
    return f"{r/10**(nd-1):.{nd-1}f}\\times10^{{{e}}}"
RC = ["R1", "R2", "R3", "R4", "R5"]; TG = [1e-3, 1e-4, 1e-5, 1e-6]; tn = {1e-3: "3", 1e-4: "4", 1e-5: "5", 1e-6: "6"}
# ---------------- (1) certified interval at d=9, p=1e-4
CS = "results/cert_d9_result.json"; C = json.load(open(f"{R}/cert_d9_result.json")); A0 = C["A0"]; AF = C["assumption_free"]
assert C["d"] == 9 and C["p"] == 1e-4 and C["plan"]["A0"] == [1, 2, 3, 4] and A0["delta"] == 0.05 and AF["delta"] == 0.05
assert sorted(map(int, A0["per"])) == list(range(5, 13)) and sorted(map(int, AF["per"])) == list(range(1, 13))
assert all(C["nfail"][str(w)] == 0 and C["nsamples"][str(w)] == 10_000_000 for w in range(1, 5))
# full-precision recomputation of the endpoints from the stored per-stratum limits (beyond-120 remainder: Chernoff bound 1e-241, below double precision, added symbolically)
_L = sum(x["P"] * x["lcb"] for x in A0["per"].values()); _U = sum(x["P"] * x["ucb"] for x in A0["per"].values()) + A0["tail"] + 1e-241
assert abs(_L - A0["L"]) <= 1e-9 * A0["L"] and abs(_U - A0["U"]) <= 1e-9 * A0["U"], (_L, A0["L"], _U, A0["U"])
_UF = sum(x["P"] * x["ucb"] for x in AF["per"].values()) + AF["tail"] + 1e-241; assert abs(_UF - AF["U"]) <= 1e-9 * AF["U"]
put("c9.L", sci_out(A0["L"], False), CS, "certified lower bound under A-0, rounded outward (down) to 3 s.f."); put("c9.U", sci_out(A0["U"], True), CS, "certified upper bound under A-0, includes the beyond-120 remainder, rounded outward (up) to 3 s.f.")
put("c9.L12", f"{math.floor(A0['L']*1e14)/100:.2f}", CS, "L in units of 1e-12, rounded down"); put("c9.U12", f"{math.ceil(A0['U']*1e13)/10:.1f}", CS, "U in units of 1e-12, rounded up"); put("c9.UoverL", f"{A0['U']/A0['L']:.0f}", CS, "U/L under A-0"); put("c9.est", sci(A0["est"]), CS, "plug-in point estimate (no guarantee)")
put("c9.free.L", sci_out(AF["L"], False), CS, "lower bound, no assumption on strata 1-4, rounded down"); put("c9.free.U", sci_out(AF["U"], True), CS, "upper bound, no assumption on strata 1-4, rounded up")
put("c9.events.A0", A0["n_events"], CS, "one-sided events in the union bound"); put("c9.events.free", AF["n_events"], CS, "")
put("c9.tail", sci(A0["tail"]), CS, "P(W>12), added to U")
nd = sum(C["nsamples"][str(w)] for w in range(5, 13)); assert C["plan"]["Ntot"] == 10 ** 9 and 0 <= nd - 10 ** 9 < 10 ** 4  # 1,000,003,205: the two smallest strata are floored at 2000 draws
put("c9.draws", sci(nd), CS, "fresh draws over strata 5-12"); put("c9.lowdraws", "10^{7}", CS, "fresh draws in each of strata 1-4"); put("c9.N", C["plan"]["N"], CS, "fault mechanisms")
for w in (5, 6, 7):
    s = A0["per"][str(w)]
    put(f"c9.w{w}.lo", sci(s["P"] * s["lcb"]), CS, "stratum contribution to L"); put(f"c9.w{w}.hi", sci(s["P"] * s["ucb"]), CS, "stratum contribution to U")
put("c9.w5.nfail", A0["per"]["5"]["nfail"], CS, "failing draws, stratum 5"); put("c9.w5.n", sci(A0["per"]["5"]["n"]), CS, "draws, stratum 5")
put("c9.w5.shareU", pct(A0["per"]["5"]["P"] * A0["per"]["5"]["ucb"] / A0["U"], 0), CS, "share of U from stratum 5, percent")
for w in (1, 2): put(f"c9.free.w{w}.hi", sci(AF["per"][str(w)]["P"] * AF["per"][str(w)]["ucb"]), CS, "contribution of a zero-failure stratum to the assumption-free U")
put("c9.free.ucb", sci(max(AF["per"][str(w)]["ucb"] for w in range(1, 5))), CS, "largest upper limit on f_w for w<=4 from 1e7 failure-free draws")
with open("tab_cert_d9.tex", "w") as f:
    f.write("\\begin{tabular}{rrrrrr}\\hline\\hline\n$w$ & $P(W=w)$ & draws & failing draws & contribution to $L$ & contribution to $U$\\\\\\hline\n")
    for w in range(5, 13):
        s = A0["per"][str(w)]; lo = s["P"] * s["lcb"]
        f.write(f"{w} & ${sci(s['P'])}$ & ${sci(s['n'])}$ & {s['nfail']} & {'$'+sci(lo)+'$' if lo > 0 else '0'} & ${sci(s['P']*s['ucb'])}$\\\\\n")
    f.write(f"$>12$ & ${sci(A0['tail'])}$ & 0 & & 0 & ${sci(A0['tail'])}$\\\\\\hline\ntotal (full precision, rounded outward) & & ${sci(nd)}$ & & ${sci_out(A0['L'], False)}$ & ${sci_out(A0['U'], True)}$\\\\\n\\hline\\hline\n\\end{{tabular}}\n".replace("{{tabular}}", "{tabular}"))
VS = "results/cert_d9_vs_recipes.json (code/cert_d9_vs_recipes.py)"; VR = json.load(open(f"{R}/cert_d9_vs_recipes.json"))
assert abs(VR["A0"]["L"] - A0["L"]) < 1e-25 and abs(VR["A0"]["U"] - A0["U"]) < 1e-25 and VR["decisive_bar_A0_R1_R5_ge_285of300"] is False
miss = {rc: VR["A0"]["recipes"][rc]["certified_missed_high"] for rc in RC}
assert miss == dict(R1=300, R2=0, R3=0, R4=239, R5=208) and all(VR["A0"]["recipes"][rc]["certified_missed_low"] == 0 and VR["A0"]["recipes"][rc]["n"] == 300 for rc in RC)
assert all(v["certified_missed_high"] == 0 and v["certified_missed_low"] == 0 for v in VR["assumption_free"]["recipes"].values())
put("c9.bar", 285, "PLAN.md (certified interval amendment)", "registered bar: replicates certified missed, of 300, for both R1 and R5")
for rc in RC:
    v = VR["A0"]["recipes"][rc]; put(f"c9.miss.{rc}", v["certified_missed_high"], VS, "replicates (of 300) whose whole 95% interval lies above U")
    f12 = lambda x: f"{x*1e12:.1f}" if x < 1e-10 else f"{x*1e12:.0f}"
    put(f"c9.int12.lo.{rc}", f12(v["median_lo"]), VS, "median lower end, units of 1e-12"); put(f"c9.int12.hi.{rc}", f12(v["median_hi"]), VS, "median upper end, units of 1e-12")
    put(f"c9.int.lo.{rc}", sci(v["median_lo"]), VS, "median lower end of the recipe's interval"); put(f"c9.int.hi.{rc}", sci(v["median_hi"]), VS, "median upper end")
put("c9.miss.R5.pct", pct(miss["R5"] / 300, 0), VS, "percent"); put("c9.miss.R4.pct", pct(miss["R4"] / 300, 0), VS, "percent")
put("c9.R1.lo.overU", f"{VR['A0']['recipes']['R1']['median_lo']/A0['U']:.1f}", VS, "median lower end of R1's interval / U")
# median predictions, recomputed exactly as make_numbers.py does for h.pred.*, and the old holdout point estimate
reps = load_reps(f"{R}/fit_d357.json"); r9 = load_reps(f"{R}/fit_d9.json"); allr = [{**a, **b} for a, b in zip(reps, r9)]; assert len(allr) == 300
hm = {rc: [] for rc in RC}
for data in allr:
    o = fit_all({d: data[d] for d in (3, 5, 9)}, 1e-4)
    for rc in RC: hm[rc].append(o[(rc, 9)][0])
PS_ = "results/fit_d357.json, fit_d9.json (code/coverage.py) and results/cert_d9_result.json"
for rc in RC:
    med = float(np.exp(np.median(hm[rc])))
    put(f"c9.pred.overU.{rc}", f"{med/A0['U']:.1f}" if med / A0["U"] < 10 else f"{med/A0['U']:.0f}", PS_, "median prediction / certified upper bound U (smallest over-prediction consistent with the interval)")
    put(f"c9.pred.overL.{rc}", f"{med/A0['L']:.1f}" if med / A0["L"] < 10 else f"{med/A0['L']:.0f}", PS_, "median prediction / certified lower bound L (largest)")
    put(f"c9.pred.inside.{rc}", sum(A0["L"] <= math.exp(x) <= A0["U"] for x in hm[rc]), PS_, "replicate point predictions inside [L,U], of 300")
H = [json.load(open(f)) for f in sorted(glob.glob(f"{R}/holdout_truth_d9_p1e-4_s*.json"))]; assert len(H) == 6
w = np.array([1 / r["se"] ** 2 for r in H]); t = np.array([r["truth"] for r in H]); ht = float((w * t).sum() / w.sum())
assert ht < A0["L"], "old holdout point estimate lies below the certified lower bound"
assert sorted(t)[4] < A0["L"] and sorted(t)[5] > A0["U"], "five runs below L, the sixth above U"
put("c9.hold.overL", f"{ht/A0['L']:.2f}", "results/holdout_truth_d9_p1e-4_s*.json + cert_d9_result.json", "old holdout point estimate / certified lower bound")
put("c9.hold.below", sum(x < A0["L"] for x in t), "same", "holdout runs (of 6) whose estimate is below the certified lower bound")
put("c9.est.overhold", f"{A0['est']/ht:.1f}", "same", "plug-in estimate of the new measurement / old holdout point estimate")
with open("tab_cert_recipes.tex", "w") as f:
    f.write("\\begin{tabular}{lcrrr}\\hline\\hline\nrecipe & interval ($10^{-12}$) & missed & pred./$U$ & pred./$L$\\\\\\hline\n")
    for rc in RC: f.write(f"{rc} & \\V{{c9.int12.lo.{rc}}}--\\V{{c9.int12.hi.{rc}}} & \\V{{c9.miss.{rc}}} & \\V{{c9.pred.overU.{rc}}} & \\V{{c9.pred.overL.{rc}}}\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
# ---------------- (2) R4' and R2-lin on the window below the crossing
RS = "results/L_r4p_bw.json (code/r4p_bw.py)"; Q = json.load(open(f"{R}/L_r4p_bw.json")); assert Q["_sanity_mismatch"] == 0
def cell(rc, d, p, lab="primary"):
    for k, v in Q.items():
        if k == "_sanity_mismatch": continue
        a, b, c, e = k.split("|")
        if a == rc and int(b) == d and abs(float(c) / p - 1) < 1e-9 and e == lab: return v
    raise KeyError((rc, d, p, lab))
for rc, key in (("R4p", "r4b"), ("R2lin", "r2b")):
    for d in (3, 5, 7):
        for p in TG:
            c = cell(rc, d, p); assert c["reps"] == 300
            put(f"{key}.cov.d{d}.t{tn[p]}", pct(c["cov"]), RS, "coverage, percent of 300 (primary)"); put(f"{key}.ratio.d{d}.t{tn[p]}", f"{c['med_ratio']:.2f}", RS, "median prediction / reference")
            if d == 7: put(f"{key}.hitmin.d7.t{tn[p]}", pct(c["hit_min"], 0), RS, "smallest hit fraction over reference +-1.96 se"); put(f"{key}.hitmax.d7.t{tn[p]}", pct(c["hit_max"], 0), RS, "largest")
low = [cell("R4p", d, p) for d in (3, 5) for p in TG]; inb = [c for c in low if c["band"] == "inband"]; ab = [c for c in low if c["band"] == "above97"]
assert len(low) == 8 and len(inb) == 5 and len(ab) == 3 and all(c["cov"] >= 0.9 for c in low)
assert all(cell("R4p", 3, p)["band"] == "above97" for p in TG[1:]) and all(cell("R4p", 5, p)["band"] == "inband" for p in TG)
put("r4b.inband", len(inb), RS, "R4' d=3,5 cells in 93-97%"); put("r4b.cells", len(low), RS, ""); put("r4b.bar", 6, "PLAN.md (amendment R4P-BW)", "registered bar")
put("r4b.above", len(ab), RS, "cells above 97%"); put("r4b.above.covmin", pct(min(c["cov"] for c in ab)), RS, ""); put("r4b.above.covmax", pct(max(c["cov"] for c in ab)), RS, "")
put("r4b.lowd.covmin", pct(min(c["cov"] for c in low)), RS, ""); put("r4b.lowd.covmax", pct(max(c["cov"] for c in low)), RS, "")
put("r4b.d5.covmin", pct(min(cell("R4p", 5, p)["cov"] for p in TG)), RS, ""); put("r4b.d5.covmax", pct(max(cell("R4p", 5, p)["cov"] for p in TG)), RS, "")
put("r4b.lowd.ratmin", f"{min(c['med_ratio'] for c in low):.2f}", RS, ""); put("r4b.lowd.ratmax", f"{max(c['med_ratio'] for c in low):.2f}", RS, "")
sl = [cell("R4p", d, p, "scaled") for d in (3, 5) for p in TG]
put("r4b.sca.inband", sum(c["band"] == "inband" for c in sl), RS, "scaled sensitivity: cells in band"); put("r4b.sca.above", sum(c["band"] == "above97" for c in sl), RS, "")
assert V["r4b.sca.inband"][0] == "4" and V["r4b.sca.above"][0] == "4"
d7 = [cell("R4p", 7, p) for p in TG]; assert all(c["band"] == "inband" and len(c["classes"]) > 1 for c in d7)
put("r4b.d7.inband", 4, RS, "d=7 cells in band at the point reference"); put("r4b.d7.undecided", sum(len(c["classes"]) > 1 for c in d7), RS, "d=7 cells whose class changes within reference +-1.96 se")
put("r4b.d7.covmin", pct(min(c["cov"] for c in d7)), RS, ""); put("r4b.d7.covmax", pct(max(c["cov"] for c in d7)), RS, "")
put("r4b.d7.ratmin", f"{min(c['med_ratio'] for c in d7):.2f}", RS, ""); put("r4b.d7.ratmax", f"{max(c['med_ratio'] for c in d7):.2f}", RS, "")
put("r4b.d7.hitmin", pct(min(c["hit_min"] for c in d7), 0), RS, ""); put("r4b.d7.relse.min", pct(min(c["truth_rel_se"] for c in d7), 0), RS, ""); put("r4b.d7.relse.max", pct(max(c["truth_rel_se"] for c in d7), 0), RS, "")
l2 = [cell("R2lin", d, p) for d in (3, 5) for p in TG]; l7 = [cell("R2lin", 7, p) for p in TG]
assert sum(c["band"] == "inband" for c in l2) == 6 and sum(c["band"] == "inband" for c in l7) == 1 and sum(c["cov"] < 0.9 for c in l7) == 3
put("r2b.inband", 6, RS, "R2-lin d=3,5 cells in band"); put("r2b.lowd.covmin", pct(min(c["cov"] for c in l2)), RS, ""); put("r2b.lowd.covmax", pct(max(c["cov"] for c in l2)), RS, "")
put("r2b.d7.inband", 1, RS, ""); put("r2b.d7.below90", 3, RS, ""); put("r2b.d7.covmin", pct(min(c["cov"] for c in l7)), RS, ""); put("r2b.d7.covmax", pct(max(c["cov"] for c in l7)), RS, "")
with open("tab_r4p_bw.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrr}\\hline\\hline\nrecipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for d in (3, 5, 7): f.write(f"R4$'$ & {d} & " + " & ".join(f"\\V{{r4b.cov.d{d}.t{t}}} (\\V{{r4b.ratio.d{d}.t{t}}})" for t in "3456") + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
with open("tab_r4p_bw_supp.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrr}\\hline\\hline\nquantity & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    f.write("R4$'$, hit fraction, smallest--largest & 7 & " + " & ".join(f"\\V{{r4b.hitmin.d7.t{t}}}--\\V{{r4b.hitmax.d7.t{t}}}" for t in "3456") + "\\\\\n")
    for d in (3, 5, 7): f.write(f"R2-lin, coverage (ratio) & {d} & " + " & ".join(f"\\V{{r2b.cov.d{d}.t{t}}} (\\V{{r2b.ratio.d{d}.t{t}}})" for t in "3456") + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
# ---------------- (3) rolling-origin discrepancy inflation (R2, R3; R1 as side check)
IS_ = "results/L_inflate_coverage.csv, results/L_inflate_verdicts.txt (code/inflate_rolling.py)"
rows = [r for r in csv.DictReader(open(f"{R}/L_inflate_coverage.csv")) if r["variant"] == "CONFIRMATORY"]
txt = open(f"{R}/L_inflate_verdicts.txt").read(); assert "BW 5/24, STRADDLE 3/24 -> FAIL" in txt
tot = {}
for win, key in (("BW", "bw"), ("STRADDLE", "st")):
    n = 0
    for rc in ("R2", "R3", "R1"):
        c = [r for r in rows if r["window"] == win and r["recipe"] == rc]; assert len(c) == 12
        cov = [float(r["cov"]) for r in c]; un = [float(r["cov_uninflated"]) for r in c]
        ib = sum(0.93 <= x <= 0.97 for x in cov); lo = sum(x < 0.9 for x in cov); hi = sum(x > 0.97 for x in cov)
        m = re.search(rf"\[{win} CONFIRMATORY\] {rc}: inflated in-band\(93-97\) (\d+)/12, <90 (\d+), >97 (\d+), median cov ([\d.]+); uninflated in-band \d+, <90 \d+, median cov ([\d.]+); median width ratio lnscale ([\d.]+)x", txt)
        assert m and (ib, lo, hi) == tuple(int(m.group(i)) for i in (1, 2, 3)), (win, rc, ib, lo, hi)
        assert abs(float(np.median(cov)) - float(m.group(4))) < 6e-4 and abs(float(np.median(un)) - float(m.group(5))) < 6e-4
        put(f"inf.{key}.{rc}.inband", ib, IS_, "cells of 12 in 93-97% after inflation"); put(f"inf.{key}.{rc}.below90", lo, IS_, ""); put(f"inf.{key}.{rc}.above97", hi, IS_, "")
        put(f"inf.{key}.{rc}.medcov", pct(float(np.median(cov)), 0), IS_, "median coverage after inflation, percent"); put(f"inf.{key}.{rc}.medcov0", pct(float(np.median(un)), 0), IS_, "before")
        put(f"inf.{key}.{rc}.width", f"{float(m.group(6)):.1f}" if float(m.group(6)) < 10 else f"{float(m.group(6)):.0f}", "results/L_inflate_verdicts.txt", "median ratio of inflated to original interval width (log scale)")
        if rc != "R1": n += ib
    tot[key] = n; put(f"inf.{key}.inband", n, IS_, "R2+R3 cells of 24 in band")
assert tot == dict(bw=5, st=3)
put("inf.bar", 18, "PLAN.md (Amendment 2)", "registered bar, cells of 24, in each window"); put("inf.cells", 24, "", "")
m = re.search(r"\[BW CONFIRMATORY\] far targets \(<=1e-4\) R2\+R3 in band: (\d+)/18, <90: (\d+)", txt); assert m and m.groups() == ("3", "11")
put("inf.bw.far.inband", 3, "results/L_inflate_verdicts.txt", ""); put("inf.bw.far.below90", 11, "results/L_inflate_verdicts.txt", ""); put("inf.bw.far.cells", 18, "", "")
with open("tab_inflate.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrrr}\\hline\\hline\nwindow & recipe & in 93--97\\% & below 90\\% & above 97\\% & median coverage, after (before) & width ratio\\\\\\hline\n")
    for key, name in (("bw", "below the crossing"), ("st", "registered")):
        for rc in ("R2", "R3", "R1"):
            f.write(f"{name} & {rc}{' (side check)' if rc == 'R1' else ''} & \\V{{inf.{key}.{rc}.inband}} & \\V{{inf.{key}.{rc}.below90}} & \\V{{inf.{key}.{rc}.above97}} & \\V{{inf.{key}.{rc}.medcov}} (\\V{{inf.{key}.{rc}.medcov0}}) & \\V{{inf.{key}.{rc}.width}}\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
with open("numbers_v5.tex", "w") as f:
    f.write("% GENERATED by make_numbers_v5.py. Do not edit. Input after numbers_v4.tex.\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers_v5.json", "w"), indent=1)
for k, (v, s, h) in V.items(): print(k, v)
