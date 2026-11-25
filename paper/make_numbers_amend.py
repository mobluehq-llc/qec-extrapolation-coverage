"""Numbers for the CORRECTIVE amendment (R4', R2-lin, R5 decomposition, reference-uncertainty bounds, single-point Lambda, d* table).
Run from this folder after make_numbers.py:  python3 make_numbers_amend.py
Reads ../results/L_amend_gemini.json (written by the registered ../code/amend_gemini.py) and ../results/fit_d357.json. Nothing is simulated or refitted
here except the finite-size crossing marker (log-linear interpolation of window medians, as in ../code/fig_amend_gemini.py).
Writes numbers_amend.tex (keys am.*), numbers_amend.json (key -> value, source) and tab_am_*.tex."""
import sys, json, math, re
import numpy as np
sys.path.insert(0, "../code")
from coverage import load_reps, PS
A = json.load(open("../results/L_amend_gemini.json")); SRC = "results/L_amend_gemini.json"
assert A["sanity_mismatch"] == 0
V = {}
def put(k, v, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), how)
def sci(x, nd=3):
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
def pct(x, nd=1): return f"{100*x:.{nd}f}"
def fx(x): return f"{x:.2f}" if x < 10 else f"{x:.1f}" if x < 100 else f"{x:.0f}"
DS = [3, 5, 7]; TG = ["0.001", "0.0001", "1e-05", "1e-06"]; tn = dict(zip(TG, "3456"))
C = A["coverage"]; RU = A["reference_uncertainty"]
for rc, tag in (("R4p", "R4p"), ("R2lin", "R2lin"), ("R4", "R4")):
    for d in DS:
        for t in TG:
            for s, st in (("primary", "pri"), ("scaled", "sca")):
                x = C[f"{rc}|{d}|{t}|{s}"]
                if rc != "R4": put(f"am.cov.{st}.{tag}.d{d}.t{tn[t]}", pct(x["cov"]), "hits/300, percent")
            x = C[f"{rc}|{d}|{t}|primary"]
            if rc != "R4": put(f"am.ratio.{tag}.d{d}.t{tn[t]}", f"{x['med_ratio']:.2f}", "median prediction / reference")
            put(f"am.hw.{tag}.d{d}.t{tn[t]}", pct(math.sqrt(x["med_width_ratio"]) - 1, 0), "median 95% half-width, percent of prediction: sqrt(upper/lower)-1")
lo = [C[f"R4p|{d}|{t}|primary"]["cov"] for d in (3, 5) for t in TG]
put("am.R4p.lowd.covmin", pct(min(lo)), "R4' primary coverage, d=3,5, min"); put("am.R4p.lowd.covmax", pct(max(lo)), "max")
r = [C[f"R4p|{d}|{t}|primary"]["med_ratio"] for d in (3, 5) for t in TG]
put("am.R4p.lowd.ratmin", f"{min(r):.2f}"); put("am.R4p.lowd.ratmax", f"{max(r):.2f}")
r = [C[f"R4p|7|{t}|primary"]["med_ratio"] for t in TG]; put("am.R4p.d7.ratmin", f"{min(r):.2f}"); put("am.R4p.d7.ratmax", f"{max(r):.2f}")
r = [C[f"R4p|{d}|{t}|primary"]["med_ratio"] for d in DS for t in TG]; put("am.R4p.all.ratmin", f"{min(r):.2f}"); put("am.R4p.all.ratmax", f"{max(r):.2f}")
h = [math.sqrt(C[f"R4p|{d}|{t}|primary"]["med_width_ratio"]) - 1 for d in DS for t in TG]
put("am.R4p.hwmin", pct(min(h), 0)); put("am.R4p.hwmax", pct(max(h), 0))
h4 = [math.sqrt(C[f"R4|{d}|{t}|primary"]["med_width_ratio"]) - 1 for d in DS for t in TG]
assert max(abs(a - b) for a, b in zip(h, h4)) < 1e-9, "R4 and R4' interval widths should be identical"
r = [C[f"R2lin|7|{t}|primary"]["med_ratio"] for t in TG]; put("am.R2lin.d7.ratmin", f"{min(r):.2f}"); put("am.R2lin.d7.ratmax", f"{max(r):.2f}")
for d in DS:
    c = [C[f"R2lin|{d}|{t}|primary"]["cov"] for t in TG]; put(f"am.R2lin.d{d}.covmin", pct(min(c))); put(f"am.R2lin.d{d}.covmax", pct(max(c)))
# R4 mechanism
M = A["R4_mechanism"]; st = [M[str(d)]["R4_over_R4p_steps"] for d in DS]
for i in range(3):
    assert max(s[i] for s in st) - min(s[i] for s in st) < 0.005
    put(f"am.mech.step{i+1}", f"{st[0][i]:.2f}", "successive-decade ratio of median R4/R4', same at every d to 2 decimals")
put("am.mech.sqrt10", f"{math.sqrt(10):.2f}")
ls = [A["ratio"][f"{d}|{t}"]["log_share"] for d in DS for t in TG[1:]]; put("am.mech.logshare.min", f"{min(ls):.2f}", "min log share at 1e-4..1e-6"); put("am.mech.logshare.max", f"{max(ls):.2f}")
# R5 decomposition
R5 = A["R5_decomposition"]
for t in TG:
    x = R5[t]
    put(f"am.r5.proj.t{tn[t]}", sci(x["ref_projection"])); put(f"am.r5.dist.t{tn[t]}", f"{x['distance_factor']:.2f}", "p_L(5)^2/p_L(3) over d=7 reference")
    put(f"am.r5.med.t{tn[t]}", sci(x["R5_median"])); put(f"am.r5.tot.t{tn[t]}", fx(x["R5_over_truth"])); put(f"am.r5.fit.t{tn[t]}", fx(x["fit_input_factor"]))
    put(f"am.r5.share.t{tn[t]}", pct(math.log(x["fit_input_factor"]) / math.log(x["R5_over_truth"]), 0), "ln(fit-input factor)/ln(total), percent")
dv = [R5[t]["distance_factor"] for t in TG]; put("am.r5.dist.min", f"{min(dv):.2f}"); put("am.r5.dist.max", f"{max(dv):.2f}")
# reference-uncertainty bounds
reg = [k for k in RU if k.endswith("|primary") and k.split("|")[0] in ("R1", "R2", "R3", "R4", "R5")]
fl = [k for k in reg if RU[k]["can_flip"]]
put("am.ru.reg.cells", len(reg)); put("am.ru.reg.flip", len(fl)); put("am.ru.reg.fixed", len(reg) - len(fl))
assert all(k.split("|")[0] in ("R2", "R3") and k.split("|")[1] == "7" for k in fl)
assert not any(RU[k]["can_flip"] for k in reg if k.split("|")[0] in ("R1", "R5", "R4"))
put("am.ru.all.cells", len(RU)); put("am.ru.all.flip", sum(v["can_flip"] for v in RU.values()))
put("am.ru.flip.hitmax.min", pct(min(RU[k]["hit_max"] for k in fl)), "smallest upper bound of hit fraction among the undecided registered cells")
put("am.ru.flip.hitmax.max", pct(max(RU[k]["hit_max"] for k in fl)))
assert all(RU[k]["hit_min"] == 0 for k in fl)
put("am.ru.flip.relse.min", pct(min(RU[k]["truth_rel_se"] for k in fl), 0)); put("am.ru.flip.relse.max", pct(max(RU[k]["truth_rel_se"] for k in fl), 0))
for rc in ("R2", "R3", "R4p"):
    for t in TG:
        x = RU[f"{rc}|7|{t}|primary"]
        put(f"am.ru.{rc}.d7.t{tn[t]}.pt", pct(x["hit_at_point"])); put(f"am.ru.{rc}.d7.t{tn[t]}.min", pct(x["hit_min"])); put(f"am.ru.{rc}.d7.t{tn[t]}.max", pct(x["hit_max"]))
        put(f"am.ru.{rc}.d7.t{tn[t]}.flag", "undecided" if x["can_flip"] else "decided")
put("am.ru.fixed.lowd.hitmax", pct(max(RU[f"{rc}|{d}|{t}|primary"]["hit_max"] for rc in ("R2", "R3") for d in (3, 5) for t in TG)), "largest hit fraction of R2,R3 at d=3,5 over the reference range")
put("am.ru.R4p.lowd.flip", sum(RU[f"R4p|{d}|{t}|primary"]["can_flip"] for d in (3, 5) for t in TG), "R4' primary cells at d=3,5 that can change class")
put("am.ru.R4p.d7.flip", sum(RU[f"R4p|7|{t}|primary"]["can_flip"] for t in TG))
# single-point Lambda
L = A["lambda_proj"]
put("am.lam.cov", pct(L["cov"])); put("am.lam.hits", L["hits"]); put("am.lam.cplo", pct(L["cp_lo"])); put("am.lam.cphi", pct(L["cp_hi"]))
put("am.lam.ratio", f"{L['med_ratio']:.2f}"); put("am.lam.q25", f"{L['q25']:.2f}"); put("am.lam.q75", f"{L['q75']:.2f}"); put("am.lam.width", f"{L['med_interval_ratio']:.2f}")
put("am.lam.ref", sci(L["truth_d7"]))
# d* table
D = A["dstar"]; names = [("truth", "reference"), ("R1", "R1"), ("R2", "R2"), ("R3", "R3"), ("R4", "R4 (registered)"), ("R4p", "R4$'$ (corrective)"), ("R2lin", "R2-lin (corrective)"), ("R5", "R5")]
for k, _ in names:
    x = D[k]; put(f"am.ds.{k}.p7", sci(x["p7"])); put(f"am.ds.{k}.lam", f"{x['lam']:.2f}"); put(f"am.ds.{k}.d", x["dstar"]); put(f"am.ds.{k}.q", f"{x['qubits']:,}".replace(",", "{,}")); put(f"am.ds.{k}.st", f"{x['spacetime_qubit_rounds']:,}".replace(",", "{,}"))
    put(f"am.ds.{k}.rel", f"{x['p7']/D['truth']['p7']:.2f}", "p_L(7) at 1e-3 over the reference")
# finite-size crossing of the d=5 and d=7 window medians (log-linear interpolation; a marker, not a threshold)
reps = load_reps("../results/fit_d357.json"); Wm = {d: np.median(np.array([r[d][0] / r[d][1] for r in reps]), 0) for d in DS}
ld = np.log(Wm[5]) - np.log(Wm[7]); cr = [math.exp(math.log(PS[i]) + ld[i] / (ld[i] - ld[i + 1]) * (math.log(PS[i + 1]) - math.log(PS[i]))) for i in range(4) if ld[i] * ld[i + 1] < 0]
assert len(cr) == 1; put("am.cross", f"{cr[0]*1e3:.1f}", "d=5/d=7 window-median crossing, units of 1e-3; results/fit_d357.json")

with open("numbers_amend.tex", "w") as f:
    f.write("% GENERATED by make_numbers_amend.py. Do not edit. Input after numbers.tex.\n")
    for k, (v, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=SRC, how=h) for k, (v, h) in V.items()}, open("numbers_amend.json", "w"), indent=1)
H = ["\\begin{tabular}{llrrrr}\\hline\\hline", "recipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline"]
def tab(fn, rows):
    open(fn, "w").write("\n".join(H + rows + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
def row(label, d, pat): return f"{label} & {d} & " + " & ".join(pat.replace("#", tn[t]) for t in TG) + "\\\\"
tab("tab_am_main.tex", [row("R4 (registered)", d, f"\\V{{cov.pri.R4.d{d}.t#}} (\\V{{ratio.R4.d{d}.t#}})") for d in DS] +
    [row("R4$'$ (corrective)", d, f"\\V{{am.cov.pri.R4p.d{d}.t#}} (\\V{{am.ratio.R4p.d{d}.t#}})") for d in DS])
tab("tab_am_supp.tex", [row("R4$'$, scaled", d, f"\\V{{am.cov.sca.R4p.d{d}.t#}}") for d in DS] +
    [row("R2-lin, primary", d, f"\\V{{am.cov.pri.R2lin.d{d}.t#}} (\\V{{am.ratio.R2lin.d{d}.t#}})") for d in DS] + [row("R2-lin, scaled", d, f"\\V{{am.cov.sca.R2lin.d{d}.t#}}") for d in DS])
tab("tab_am_ru.tex", [row(lab, 7, f"\\V{{am.ru.{rc}.d7.t#.min}}--\\V{{am.ru.{rc}.d7.t#.max}}") for rc, lab in (("R2", "R2"), ("R3", "R3"), ("R4p", "R4$'$ (corrective)"))])
L = ["\\begin{tabular}{lrrrrr}\\hline\\hline", "target $p$ & projection & distance & R5 median & total & fit-input\\\\\\hline"]
for t in TG: L.append(f"$10^{{-{tn[t]}}}$ & $\\V{{am.r5.proj.t{tn[t]}}}$ & \\V{{am.r5.dist.t{tn[t]}}} & $\\V{{am.r5.med.t{tn[t]}}}$ & \\V{{am.r5.tot.t{tn[t]}}} & \\V{{am.r5.fit.t{tn[t]}}}\\\\")
open("tab_am_r5.tex", "w").write("\n".join(L + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
L = ["\\begin{tabular}{lrrrrr}\\hline\\hline", "source & $\\pL(7)$/ref. & $\\Lambda$ & $d^*$ & qubits & qubit-rounds\\\\\\hline"]
for k, lab in names: L.append(f"{lab} & \\V{{am.ds.{k}.rel}} & \\V{{am.ds.{k}.lam}} & \\V{{am.ds.{k}.d}} & \\V{{am.ds.{k}.q}} & \\V{{am.ds.{k}.st}}\\\\")
open("tab_am_dstar.tex", "w").write("\n".join(L + ["\\hline\\hline", "\\end{tabular}"]) + "\n")
print(len(V), "amendment numbers written")
for k in ["am.R4p.lowd.covmin", "am.R4p.lowd.covmax", "am.R4p.hwmin", "am.R4p.hwmax", "am.mech.step1", "am.mech.step2", "am.mech.step3", "am.r5.dist.t6", "am.r5.fit.t6", "am.r5.tot.t6", "am.r5.share.t6", "am.ru.reg.flip", "am.ru.all.flip", "am.ru.flip.hitmax.min", "am.ru.flip.hitmax.max",
          "am.lam.cov", "am.cross", "am.ds.R4.rel", "am.ds.R1.rel", "am.ds.R5.rel", "am.ds.R2.rel", "am.ru.R4p.lowd.flip", "am.ru.R4p.d7.flip", "am.mech.logshare.min", "am.R2lin.d7.ratmin", "am.R2lin.d7.ratmax", "am.R2lin.d3.covmin","am.R2lin.d3.covmax","am.R2lin.d5.covmin","am.R2lin.d5.covmax","am.R2lin.d7.covmin","am.R2lin.d7.covmax","am.R4p.d7.ratmin","am.R4p.d7.ratmax","am.R4p.lowd.ratmin","am.R4p.lowd.ratmax","am.ru.flip.relse.min","am.ru.flip.relse.max","am.hw.R4p.d3.t3","am.hw.R4p.d3.t6","am.hw.R4p.d7.t6"]: print(k, V[k][0])
