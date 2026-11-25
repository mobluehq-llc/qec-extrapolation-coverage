"""Numbers for draft v4: the below-threshold amendment (PLAN.md 'BELOW-THRESHOLD AMENDMENT', results in ../RESULTS.md, commit f52ae75)
and the corrected count of d=7 reference-estimator runs by the way they ended.
Run from this folder after make_numbers.py, make_numbers_amend.py and make_numbers_v3.py:
  python3 make_numbers_v4.py
Reads only stored files (../results/fit_bw.json, L_bw_coverage*.csv, L_coverage.csv, truth files). Nothing is simulated. The frozen
code/coverage_bw.py is imported unchanged; the coverage CSVs written on a compute host are re-derived here and asserted equal (hits per cell)
before any number is emitted. Writes numbers_v4.tex (keys bw.*, m4.*), numbers_v4.json, tab_bw_*.tex, figs/fig4_windows.{pdf,png}."""
import sys, json, csv, glob, math, re, os
import numpy as np
sys.path.insert(0, "../code")
os.chdir("../code")
from coverage_bw import fit_all, load_reps, wls, coverage, PS, X0, TARGETS
from analyze_bw import load_truth
T = load_truth(); os.chdir("../paper")
R = "../results"; V = {}
def put(k, v, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), src, how)
def pct(x, nd=1): return f"{100*x:.{nd}f}"
def fx(x): return f"{x:.2f}" if x < 10 else f"{x:.0f}"   # ratios near 1 need two decimals here
def sci(x, nd=3):
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
RC = ["R1", "R2", "R3", "R4", "R5"]; DS = [3, 5, 7]; TG = [1e-3, 1e-4, 1e-5, 1e-6]; FAR = TG[1:]
tn = {1e-3: "3", 1e-4: "4", 1e-5: "5", 1e-6: "6"}
assert list(PS) == [1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3] and TARGETS == TG
truth = {(d, p): T[(d, p)]["truth"] for d in DS for p in TG}
reps = load_reps(f"{R}/fit_bw.json"); assert len(reps) == 300
SRC = "results/fit_bw.json, truth files; code/coverage_bw.py (frozen), re-derived = results/L_bw_coverage.csv"
put("bw.reps", len(reps), "results/fit_bw.json", "replicates")
# ---- coverage, re-derived and checked against the files written on a compute host
COV = {}
for lab, scaled, fn in (("pri", False, "L_bw_coverage.csv"), ("sca", True, "L_bw_coverage_scaled.csv")):
    rows = coverage(reps, truth, scaled=scaled)
    disk = {(r["recipe"], int(r["d"]), float(r["target"])): int(r["hits"]) for r in csv.DictReader(open(f"{R}/{fn}"))}
    mine = {(r["recipe"], r["d"], r["target"]): r["hits"] for r in rows}
    assert mine == disk and len(mine) == 52, f"re-derived coverage differs from {fn}"
    COV[lab] = {(r["recipe"], r["d"], r["target"]): r for r in rows}
    for k, r in COV[lab].items(): put(f"bw.cov.{lab}.{k[0]}.d{k[1]}.t{tn[k[2]]}", pct(r["cov"]), SRC if lab == "pri" else SRC.replace("coverage.csv", "coverage_scaled.csv"), "coverage, percent of 300")
P = COV["pri"]; S = COV["sca"]
def cells(C, rc, tg=TG): return [C[k] for k in C if k[0] == rc and k[2] in tg]
# ---- registered predictions BW-1..BW-4 (counts as in results/L_bw_verdicts.txt and code/bw_extra.py)
VS = "results/L_bw_coverage.csv; rule in PLAN.md (below-threshold amendment); cf. results/L_bw_verdicts.txt"
for lab, C in (("pri", P), ("sca", S)):
    for rc in ("R1", "R5"):
        c = cells(C, rc, FAR)
        put(f"bw1.{lab}.{rc}.mis", sum(r["miscal"] for r in c), VS, "far cells (targets <= 1e-4) with coverage < 90%"); put(f"bw1.{lab}.{rc}.firm", sum(r["miscal_firm"] for r in c), VS, "of which CP upper limit < 90%")
        if lab == "pri": put(f"bw1.{rc}.cells", len(c), VS, "far cells")
assert V["bw1.pri.R1.mis"][0] == "0" and V["bw1.pri.R5.mis"][0] == "2"
vt = open(f"{R}/L_bw_verdicts.txt").read()
assert "primary(absolute sigma) H1 (R1 and R5 miscalibrated in >=2/3 of cells at >=2 decades): NOT SUPPORTED" in vt and "sensitivity(chi2/dof-scaled) H1 (R1 and R5 miscalibrated in >=2/3 of cells at >=2 decades): NOT SUPPORTED" in vt
old = {(r["recipe"], int(r["d"]), float(r["target"])): float(r["cov"]) for r in csv.DictReader(open(f"{R}/L_coverage.csv"))}
put("bw2.higher", sum(P[k]["cov"] > old[k] for k in P if k[0] == "R1"), VS + " vs results/L_coverage.csv", "R1 cells with coverage strictly above the registered-window table"); put("bw2.cells", 12, VS, "")
r15 = cells(P, "R1") + cells(P, "R5")
put("bw3.inband", sum(r["calibrated"] for r in r15), VS, "R1 and R5 cells in 93-97%"); put("bw3.cells", len(r15), VS, "")
put("bw3.inband.R1", sum(r["calibrated"] for r in cells(P, "R1")), VS, "R1 cells in 93-97%"); put("bw3.inband.R5", sum(r["calibrated"] for r in cells(P, "R5")), VS, "")
r1 = [r["cov"] for r in cells(P, "R1")]
put("bw.R1.median", pct(float(np.median(r1))), VS, "median R1 coverage over 12 cells"); put("bw.R1.covmin", pct(min(r1)), VS, ""); put("bw.R1.covmax", pct(max(r1)), VS, "")
put("bw.R1.mis", sum(r["miscal"] for r in cells(P, "R1")), VS, "R1 cells below 90%"); put("bw.R1.cplo.min", pct(min(r["cp_lo"] for r in cells(P, "R1"))), VS, "smallest CP lower limit, R1"); put("bw.R1.cphi.max", pct(max(r["cp_hi"] for r in cells(P, "R1"))), VS, "largest CP upper limit, R1")
put("bw.R1.cphi.below95", sum(r["cp_hi"] < 0.95 for r in cells(P, "R1")), VS, "R1 cells whose CP upper limit is below 95%")
r1s = [r["cov"] for r in cells(S, "R1")]
put("bw.sca.R1.inband", sum(r["calibrated"] for r in cells(S, "R1")), VS, "scaled: R1 cells in 93-97%"); put("bw.sca.R1.median", pct(float(np.median(r1s))), VS, ""); put("bw.sca.R1.covmin", pct(min(r1s)), VS, ""); put("bw.sca.R1.covmax", pct(max(r1s)), VS, "")
put("bw.sca.R5.median", pct(float(np.median([r["cov"] for r in cells(S, "R5")]))), VS, "")
for rc in ("R2", "R3"):
    put(f"bw4.{rc}.mis", sum(r["miscal"] for r in cells(P, rc)), VS, "cells below 90%"); put(f"bw4.{rc}.firm", sum(r["miscal_firm"] for r in cells(P, rc)), VS, "firmly"); put(f"bw4.{rc}.inband", sum(r["calibrated"] for r in cells(P, rc)), VS, "cells in 93-97%")
    put(f"bw4.{rc}.zero", sum(r["hits"] == 0 for r in cells(P, rc)), VS, "cells with no hit")
    put(f"bw4.sca.{rc}.mis", sum(r["miscal"] for r in cells(S, rc)), VS, "scaled: cells below 90%")
f23 = cells(P, "R2") + cells(P, "R3"); far23 = cells(P, "R2", FAR) + cells(P, "R3", FAR)
put("bw4.mis", sum(r["miscal"] for r in f23), VS, "R2+R3 cells below 90%"); put("bw4.cells", len(f23), VS, ""); put("bw4.firm", sum(r["miscal_firm"] for r in f23), VS, "")
put("bw4.far.mis", sum(r["miscal"] for r in far23), VS, ""); put("bw4.far.cells", len(far23), VS, ""); put("bw4.notmis", len(f23) - sum(r["miscal"] for r in f23), VS, "R2+R3 cells at or above 90%")
put("bw4.sca.mis", int(V["bw4.sca.R2.mis"][0]) + int(V["bw4.sca.R3.mis"][0]), VS, "scaled: R2+R3 cells below 90%")
assert (V["bw2.higher"][0], V["bw3.inband"][0], V["bw.R1.median"][0], V["bw4.mis"][0], V["bw4.far.mis"][0]) == ("12", "5", "92.5", "19", "16")
r4 = cells(P, "R4"); put("bw.R4.t3.covmin", pct(min(r["cov"] for r in r4 if r["target"] == 1e-3)), VS, ""); put("bw.R4.t3.covmax", pct(max(r["cov"] for r in r4 if r["target"] == 1e-3)), VS, "")
put("bw.R4.far.covmax", pct(max(r["cov"] for r in r4 if r["target"] < 1e-3)), VS, "largest R4 coverage at targets <= 1e-4"); put("bw.R4.mis", sum(r["miscal"] for r in r4), VS, "")
r5 = cells(P, "R5"); put("bw.R5.covmin", pct(min(r["cov"] for r in r5)), VS, ""); put("bw.R5.covmax", pct(max(r["cov"] for r in r5)), VS, ""); put("bw.R5.mis.all", sum(r["miscal"] for r in r5), VS, "R5 cells below 90%, all four targets")
allc = list(P.values()); put("bw.cells", len(allc), VS, ""); put("bw.cells.mis", sum(r["miscal"] for r in allc), VS, "cells below 90%"); put("bw.cells.zero", sum(r["hits"] == 0 for r in allc), VS, "cells with no hit"); put("bw.cells.inband", sum(r["calibrated"] for r in allc), VS, "")
# ---- direction and interval width (descriptive; same fits)
LR = {}; HW = {}
for data in reps:
    for p in TG:
        for (rc, d), (mu, var) in fit_all({d: data[d] for d in DS}, p).items():
            LR.setdefault((rc, d, p), []).append(mu - math.log(truth[(d, p)])); HW.setdefault((rc, d, p), []).append(1.96 * math.sqrt(var))
RAT = {k: float(np.exp(np.median(v))) for k, v in LR.items()}
RS = "results/fit_bw.json, truth files; code/coverage_bw.py fit_all (as code/bw_extra.py)"
for (rc, d, p), v in RAT.items(): put(f"bw.ratio.{rc}.d{d}.t{tn[p]}", fx(v), RS, "median over 300 replicates of prediction/reference")
for rc in RC:
    v = [RAT[k] for k in RAT if k[0] == rc]; put(f"bw.ratio.{rc}.min", fx(min(v)), RS, "min over cells"); put(f"bw.ratio.{rc}.max", fx(max(v)), RS, "max over cells")
    hw = [float(np.median(HW[k])) for k in HW if k[0] == rc]
    put(f"bw.hw.{rc}.min", pct(math.exp(min(hw)) - 1, 0 if math.exp(min(hw)) - 1 > 0.2 else 1), RS, "smallest cell-median 95% half-width, percent of prediction"); put(f"bw.hw.{rc}.max", pct(math.exp(max(hw)) - 1, 0 if math.exp(max(hw)) - 1 > 0.2 else 1), RS, "largest")
for p in TG:
    hw = [float(np.median(HW[("R1", d, p)])) for d in DS]; put(f"bw.hw.R1.t{tn[p]}.min", pct(math.exp(min(hw)) - 1, 0), RS, "R1 cell-median half-width at this target, min over d, percent"); put(f"bw.hw.R1.t{tn[p]}.max", pct(math.exp(max(hw)) - 1, 0), RS, "max over d")
v = [RAT[(rc, d, p)] for rc in ("R2", "R3") for d in DS for p in TG]; put("bw.ratio.fixed.min", fx(min(v)), RS, ""); put("bw.ratio.fixed.max", fx(max(v)), RS, "")
put("bw.fixed.low", sum(x < 1 for x in v), RS, "R2,R3 cells with median ratio below 1"); put("bw.fixed.maxdev", f"{100*max(abs(x-1) for x in v):.0f}", RS, "largest |median ratio - 1|, percent")
for t in "46": put(f"bw.ratio.R4.t{t}.min", fx(min(RAT[("R4", d, 10.0 ** -int(t))] for d in DS)), RS, "over d"); put(f"bw.ratio.R4.t{t}.max", fx(max(RAT[("R4", d, 10.0 ** -int(t))] for d in DS)), RS, "over d")
for i, (a, b) in enumerate(((1e-4, 1e-5), (1e-5, 1e-6)), 1):
    s = [RAT[("R4", d, b)] / RAT[("R4", d, a)] for d in DS]; put(f"bw.R4.step{i}.min", f"{min(s):.1f}", RS, "growth of the R4 median ratio over one decade, min over d"); put(f"bw.R4.step{i}.max", f"{max(s):.1f}", RS, "max over d")
put("bw.ratio.R1.lowd.t6.min", fx(min(RAT[("R1", d, 1e-6)] for d in (3, 5))), RS, ""); put("bw.ratio.R1.lowd.t6.max", fx(max(RAT[("R1", d, 1e-6)] for d in (3, 5))), RS, "")
NN = json.load(open("numbers.json")); put("bw.old.ratio.R1.d7.t6", NN["ratio.R1.d7.t6"]["value"], "numbers.json (registered window)", "R1 median ratio, d=7, 1e-6, registered window")
assert abs(RAT[("R1", 7, 1e-6)] - 0.85) < 0.005
for d in DS:
    sl = []; c2 = []
    for r in reps:
        k, n = r[d]; Pq = k / n; b, C, chi2, dof = wls(np.c_[np.ones(5), X0], np.log(Pq), k / (1 - Pq)); sl.append(b[1]); c2.append(chi2 / dof)
    put(f"bw.slope.d{d}", f"{np.median(sl):.2f}", RS, "median fitted R1 exponent, below-crossing window"); put(f"bw.chidof.d{d}", f"{np.median(c2):.1f}", RS, "median chi2/dof of the R1 fit")
    Pm = np.array([r[d][0] / r[d][1] for r in reps]).mean(0); nm = np.array([r[d][1] for r in reps]).mean(0)
    for ip in range(5): put(f"bw.win.d{d}.p{ip}", sci(Pm[ip]), "results/fit_bw.json", "mean over replicates of k/n")
    put(f"bw.shots.d{d}.p0", sci(nm[0], 2), "results/fit_bw.json", "mean shots per replicate at the lowest window point")
    assert all((r[d][0] >= 300).all() for r in reps)
put("bw.slope.maxdev", f"{max(abs(float(V[f'bw.slope.d{d}'][0])-(d+1)/2) for d in DS):.2f}", RS, "largest |median R1 exponent - (d+1)/2|")
# ---- window geometry (arithmetic)
put("bw.px", "8.28", "PLAN.md (below-threshold amendment): crossing of d=3 and d=7 window means, results/fit_d357_1121.json", "x 1e-3"); put("bw.topfrac", f"{4e-3/8.28e-3:.2f}", "arithmetic", "top window point / p_x")
for lab, lo in (("new", 1.5e-3), ("old", 4e-3)):
    put(f"bw.gap.{lab}.t3", f"{math.log10(lo/1e-3):.1f}", "arithmetic", "decades from lower window edge to target 1e-3"); put(f"bw.gap.{lab}.t4", f"{math.log10(lo/1e-4):.1f}", "arithmetic", "to 1e-4"); put(f"bw.gap.{lab}.t6", f"{math.log10(lo/1e-6):.1f}", "arithmetic", "to 1e-6")
put("bw.gapdiff", f"{math.log10(4e-3/1.5e-3):.1f}", "arithmetic", "decades by which the lower edge of the second window is nearer the targets")
put("bw.span.new", f"{4e-3/1.5e-3:.1f}", "arithmetic", "ratio of highest to lowest window point"); put("bw.span.old", f"{1.5e-2/4e-3:.1f}", "arithmetic", "")
# ---- d=7 reference-estimator runs, by the way they ended (corrects v3 wording 'every run at d=7 reached its time limit')
runs = [json.load(open(f)) for f in sorted(glob.glob(f"{R}/truth_d7_p*.json")) + sorted(glob.glob(f"{R}/truthextra_d7_p*.json"))]
def ended(x): return "rule" if x["se"] <= 0.02 * x["truth"] else ("draws" if x["totM"] > 40_000_000 and x["secs"] <= 2400 else ("time" if x["secs"] > 2400 else "?"))
E = [(x["p"], ended(x)) for x in runs]; assert all(e != "?" for _, e in E)
MS = "results/truth_d7_p*.json + truthextra_d7_p*.json (fields se, truth, secs, totM); stop order in code/truth.py; 'STOP by budget' lines in logs/L_cov_manifest.jsonl (4 second runs) and in the first-run logs on a compute host (read 2026-11-24)"
put("m4.d7.runs", len(runs), MS, "all stored d=7 runs of the reference estimator"); put("m4.d7.time", sum(e == "time" for _, e in E), MS, "ended at the 2400 s limit"); put("m4.d7.draws", sum(e == "draws" for _, e in E), MS, "ended at the 4e7-draw limit before 2400 s"); put("m4.d7.rule", sum(e == "rule" for _, e in E), MS, "stopped by the 2% rule")
ref = [e for p, e in E if any(abs(p / t - 1) < 1e-9 for t in TG)]; put("m4.d7.ref.runs", len(ref), MS, "runs behind the four d=7 reference cells"); put("m4.d7.ref.time", sum(e == "time" for e in ref), MS, "of which ended at the time limit")
put("m4.d7.val.runs", len(runs) - len(ref), MS, "d=7 runs used only in the check against direct Monte Carlo")
assert [V[k][0] for k in ("m4.d7.runs", "m4.d7.time", "m4.d7.draws", "m4.d7.rule", "m4.d7.ref.runs", "m4.d7.ref.time")] == ["12", "9", "1", "2", "8", "8"]
dr = [x for x in runs if ended(x) == "draws"][0]; put("m4.d7.draws.p", f"{dr['p']*1e3:.1f}", MS, "p of the run that ended at the draw limit, x 1e-3"); put("m4.d7.draws.relse", pct(dr["se"] / dr["truth"]), MS, "its se/estimate"); put("m4.d7.draws.secs", f"{dr['secs']:.0f}", MS, "")
tv = [x for x in runs if ended(x) == "time" and not any(abs(x["p"] / t - 1) < 1e-9 for t in TG)]; assert len(tv) == 1
put("m4.d7.timeval.p", f"{tv[0]['p']*1e3:.0f}", MS, "p of the validation-only run that ended at the time limit, x 1e-3"); put("m4.d7.timeval.relse", pct(tv[0]["se"] / tv[0]["truth"]), MS, "")
# (stale guard removed for the bundle: this script predates the d=9 certified result, which its outputs do not use)
# ---- tables
hdr = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
def rows(fmt): return [f"{rc} & {d} & " + " & ".join(fmt(rc, d, tn[p]) for p in TG) + "\\\\" for rc in RC for d in DS if (rc, d, 1e-3) in P]
open("tab_bw_main.tex", "w").write("\n".join(hdr + rows(lambda rc, d, t: f"\\V{{bw.cov.pri.{rc}.d{d}.t{t}}} (\\V{{bw.ratio.{rc}.d{d}.t{t}}})") + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
open("tab_bw_scaled.tex", "w").write("\n".join(hdr + rows(lambda rc, d, t: f"\\V{{bw.cov.sca.{rc}.d{d}.t{t}}}") + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
L = ["\\begin{tabular}{lrrrrr}\\hline\\hline", "$d$ & $1.5\\times10^{-3}$ & $2\\times10^{-3}$ & $2.5\\times10^{-3}$ & $3\\times10^{-3}$ & $4\\times10^{-3}$\\\\\\hline"]
L += [f"{d} & " + " & ".join(f"$\\V{{bw.win.d{d}.p{ip}}}$" for ip in range(5)) + "\\\\" for d in DS]
open("tab_bw_window.tex", "w").write("\n".join(L + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
with open("numbers_v4.tex", "w") as f:
    f.write("% GENERATED by make_numbers_v4.py. Do not edit. Input after numbers_v3.tex.\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers_v4.json", "w"), indent=1)
# ---- figure (plotting only): R1 on the two windows
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit(0)   # figure only; every number and table is already written above
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
oldrat = {(d, p): float(NN[f"ratio.R1.d{d}.t{tn[p]}"]["value"]) for d in DS for p in TG}
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7)); MK = {3: "o", 5: "s", 7: "^"}
for d in DS:
    ax[0].plot(TG, [100 * old[("R1", d, p)] for p in TG], MK[d] + "--", color="#0072B2", mfc="none", ms=4.5, lw=0.8, label=f"registered window, $d={d}$")
    ax[0].plot(TG, [100 * P[("R1", d, p)]["cov"] for p in TG], MK[d] + "-", color="#D55E00", ms=4.5, lw=0.8, label=f"window below the crossing, $d={d}$")
    ax[1].plot(TG, [oldrat[(d, p)] for p in TG], MK[d] + "--", color="#0072B2", mfc="none", ms=4.5, lw=0.8); ax[1].plot(TG, [RAT[("R1", d, p)] for p in TG], MK[d] + "-", color="#D55E00", ms=4.5, lw=0.8)
ax[0].axhspan(93, 97, color="0.85", lw=0); ax[0].axhline(90, color="0.4", lw=0.6, ls=":"); ax[0].set_ylim(-3, 103); ax[0].set_ylabel("coverage of R1 95% interval (%)"); h, l = ax[0].get_legend_handles_labels(); ax[1].legend(h, l, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.97))
ax[1].axhline(1, color="0.4", lw=0.6, ls=":"); ax[1].set_yscale("log"); ax[1].set_ylabel("R1 median prediction / reference")
for a, t in zip(ax, "ab"): a.set_xscale("log"); a.invert_xaxis(); a.set_xlabel("target $p$"); a.text(0.02, 1.02, f"({t})", transform=a.transAxes)
fig.tight_layout(); fig.savefig("figs/fig4_windows.pdf"); fig.savefig("figs/fig4_windows.png", dpi=200)
for k, (v, s, h) in V.items(): print(k, v)
