"""Numbers for draft v8 (final checks of v7 by Gemini and ChatGPT, 2026-11-24/25). Run from this folder (stdlib only): python3 make_numbers_v7.py
Reads results/L_final_checks_v7.json (code/final_checks_v7.py) and results/L_betting_monotone_v7.json (code/betting_monotone_check.py).
Writes numbers_v7.tex (keys v7.*), numbers_v7.json, tab_v7_cov.tex (Table III, one decimal throughout), tab_v7_cp.tex (hits and exact intervals).
One rule for robust/undecided on integer hit counts: robust-below kmax<270, at-or-above kmin>=270, else undecided."""
import json, re
R = "../results"; V = {}
def put(k, v, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), src, how)
F = json.load(open(f"{R}/L_final_checks_v7.json")); BM = json.load(open(f"{R}/L_betting_monotone_v7.json"))
RS = "results/L_final_checks_v7.json"; N = 300
TG = ["0.001", "0.0001", "1e-05", "1e-06"]; tn = dict(zip(TG, "3456"))
for key, s in F["summary"].items():
    w, rc, sub = key.split("|")
    for l in ("below", "undecided", "above"): put(f"v7.{w}.{rc}.{sub}.{l}", s[l], RS, f"{rc} cells {sub} with label {l}")
    put(f"v7.{w}.{rc}.{sub}.pt", s["point_below90"], RS, "cells below 90% at the reference point estimate"); put(f"v7.{w}.{rc}.{sub}.n", s["n"], RS, "cells")
for w in ("w1", "w2"):
    for l in ("below", "undecided", "above"):
        put(f"v7.{w}.fx.{l}", sum(F["summary"][f"{w}|{rc}|all"][l] for rc in ("R2", "R3")), RS, f"R2+R3 cells with label {l}")
cells = F["cells"]
def pc(x): return f"{100*x:.1f}"
for k, c in cells.items():
    w, rc, d, t = k.split("|")
    put(f"v7.k.{w}.{rc}.d{d}.t{tn[t]}", c["k"], RS, "hits of 300")
    put(f"v7.cv.{w}.{rc}.d{d}.t{tn[t]}", pc(c["k"]/N), RS, "coverage at the reference point estimate, percent, one decimal")
    put(f"v7.rg.{w}.{rc}.d{d}.t{tn[t]}", f"{pc(c['kmin']/N)}--{pc(c['kmax']/N)}", RS, "hit-fraction range, one decimal, from integer counts")
    put(f"v7.cp.{w}.{rc}.d{d}.t{tn[t]}", f"{pc(c['cp'][0])}--{pc(c['cp'][1])}", RS, "Clopper-Pearson 95% interval, percent")
    assert c["kmin"] <= c["k"] <= c["kmax"], k      # point inside its own range, on integers
NAME = dict(R1="R1", R2="R2", R3="R3", R4p="R4$'$", R5="R5")
ROWS = [(rc, d) for rc in ("R2", "R3", "R1", "R4p") for d in (3, 5, 7)] + [("R5", 7)]
with open("tab_v7_cov.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrrrrrr}\\hline\\hline\n & & \\multicolumn{4}{c}{first window ($4$--$15\\times10^{-3}$)} & \\multicolumn{4}{c}{second window ($1.5$--$4\\times10^{-3}$)}\\\\\nrecipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for rc, d in ROWS:
        f.write(f"{NAME[rc]} & {d}" + "".join(f" & \\V{{v7.cv.{w}.{rc}.d{d}.t{tn[t]}}} (\\V{{v7.rg.{w}.{rc}.d{d}.t{tn[t]}}})" for w in ("w1", "w2") for t in TG) + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
with open("tab_v7_cp.tex", "w") as f:
    f.write("\\begin{tabular}{llcccc}\\hline\\hline\nrecipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for w, wn in (("w1", "First window"), ("w2", "Second window")):
        f.write(f"\\multicolumn{{6}}{{l}}{{\\emph{{{wn}}}}}\\\\\n")
        for rc, d in ROWS:
            f.write(f"{NAME[rc]} & {d}" + "".join(f" & \\V{{v7.k.{w}.{rc}.d{d}.t{tn[t]}}}\\;[\\V{{v7.cp.{w}.{rc}.d{d}.t{tn[t]}}}]" for t in TG) + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
S = BM["strata"]
put("v7.w2.R3.d3.below", F["summary"]["w2|R3|lowd"]["below"] - sum(1 for t in TG if cells[f"w2|R3|5|{t}"]["label"] == "below"), RS, "R3 d=3 second-window cells robustly below 90% (derived: lowd below minus d=5 below); text says all d=3,5 shortfalls are at d=5")
assert V["v7.w2.R3.d3.below"][0] == "0"
put("v7.bet.mono", BM["worst_monotonicity_violation"], "results/L_betting_monotone_v7.json", "largest violation of monotonicity of log K+ (non-increasing) or log K- (non-decreasing) on a 4000-point grid, 8 strata")
put("v7.bet.adm", BM["min_admissible_at_bounds"], "results/L_betting_monotone_v7.json", "smallest number of admissible lower-test multipliers at any computed bound")
put("v7.bet.cut", f"{min(x['first_cutoff_m'] for x in S.values()):.2f}", "results/L_betting_monotone_v7.json", "smallest trial mean at which a multiplier becomes inadmissible, b/30")
put("v7.bet.ratio", f"{min(x['cutoff_over_max_bound'] for x in S.values()):.1f}", "results/L_betting_monotone_v7.json", "first cutoff over the largest computed bound, smallest over strata")
put("v7.bet.bf", f"{100*max(x['rel_diff_ucb'] for x in S.values()):.1f}", "results/L_betting_monotone_v7.json", "largest relative difference between the bisection upper bound and a dense-grid brute-force root, percent (grid step <= 0.9%)")
put("v7.bet.bfl", f"{100*max(x['rel_diff_lcb'] for x in S.values() if x['rel_diff_lcb'] is not None):.1f}", "results/L_betting_monotone_v7.json", "same for the lower bound where it is nonzero, percent")
sc = F["interval_scale"]
put("v7.sc.lower", f"{100*sc['lower_rel']:.0f}", RS, "lower end below the prediction at the median 219% upper half-width, percent"); put("v7.sc.ul", f"{sc['U_over_L']:.1f}", RS, "upper/lower at that median")
put("v7.sc.sym", f"{100*sc['symmetric_halfwidth_UminusL_over_2p']:.0f}", RS, "(U-L)/(2p) at that median, percent")

# ---- v9 (arXiv v2): controls A/B outcomes and the zero-fault result, parsed from stored result files
CA = open(f"{R}/L_control_a.txt").read(); CB = open(f"{R}/L_control_b.txt").read()
_do = [(float(a), int(b), int(c)) for a, b, c in re.findall(r"CONTROL-A R[1-4] dataonly: reject rate ([0-9.]+); cells with reject rate>=0.95: (\d+)/(\d+)", CA)]
assert len(_do) == 4 and all(r == 1.0 and b == c == 15 for r, b, c in _do), _do
put("ca.do.cells", 15, "results/L_control_a.txt", "cells (3 distances x 5 window points) in which every recipe R1-R4 rejects the data-only data at rate >=0.95")
put("ca.do.rate", "1.000", "results/L_control_a.txt", "reject rate of the data-only dataset, R1-R4")
_cb = re.findall(r"CONTROL-B (R[1-5]): miscalibrated (\d+)/(\d+) cells; median coverage ([0-9.]+)", CB)
assert len(_cb) == 5 and all(a == b for _, a, b, _ in _cb), _cb
assert [c for _, _, c, _ in _cb] == ["12", "12", "12", "12", "4"]
put("cb.cells.R14", 12, "results/L_control_b.txt", "cells per recipe R1-R4, all miscalibrated"); put("cb.cells.R5", 4, "results/L_control_b.txt", "cells for R5, all miscalibrated")
put("cb.maxmed", f"{100*max(float(m) for *_, m in _cb):.1f}", "results/L_control_b.txt", "largest median coverage over recipes, percent")
ZF = json.load(open(f"{R}/L_zerofault.json"))["verdicts"]
assert ZF["P1"]["MET"] is False and ZF["P2"]["MET"] is False and ZF["P3"]["MET"] is False
put("zf.p1.in", ZF["P1"]["in_band"], "results/L_zerofault.json", "cells with median prediction/reference in 0.90-1.10"); put("zf.p1.n", ZF["P1"]["cells"], "results/L_zerofault.json", "")
put("zf.p2.below", ZF["P2"]["below90"], "results/L_zerofault.json", "cells below 90% coverage"); put("zf.p3.rb", ZF["P3"]["robust_below"], "results/L_zerofault.json", "second-window R2E cells robust-below"); put("zf.p3.n", ZF["P3"]["cells"], "results/L_zerofault.json", "")
with open("numbers_v7.tex", "w") as f:
    f.write("% GENERATED by make_numbers_v7.py. Do not edit. Input after numbers_v6.tex.\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers_v7.json", "w"), indent=1)
print(len(V), "keys")
