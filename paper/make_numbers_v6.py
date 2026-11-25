"""Numbers for draft v6 (response to the Gemini and ChatGPT red-team reviews of v5; PLAN.md 'CORRECTIVE AMENDMENT (Gemini
red-team v5 checks)' and its addendum; checks committed in 71dab71, summary in ../REDTEAM_CHECKS_GEMINI.md).
Run from this folder (standard library only, no numpy):  python3 make_numbers_v6.py
Reads only stored result files; nothing is simulated or refitted and no stored result is changed:
  results/L_redteam_gemini.json  (code/redteam_checks_gemini.py): per-round normalisation, d=7 reference checks, surrogate bracket, tail, crossings
  results/L_amend_gemini.json    (code/amend_gemini.py): hit-fraction range under reference uncertainty, registered window
  results/L_r4p_bw.json          (code/r4p_bw.py): the same range on the window below the crossing, all recipes
  results/L_coverage.csv, results/L_bw_coverage.csv: frozen hit counts (asserted equal)
Writes numbers_v6.tex (keys v6.*), numbers_v6.json, tab_v6_cov.tex, tab_v6_ratio.tex, tab_v6_r2d7.tex, tab_v6_r5.tex,
tab_v6_dstar.tex, tab_v6_bracket.tex, tab_v6_d7runs.tex. Everything here is CORRECTIVE or descriptive; no frozen verdict changes."""
import json, csv, math, re
R = "../results"; V = {}
def put(k, v, src, how=""):
    assert re.fullmatch(r"[A-Za-z0-9.]+", k) and k not in V, k
    V[k] = (str(v), src, how)
def sci(x, nd=2):
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
def comma(n): return f"{n:,}".replace(",", "{,}")
TG = ["0.001", "0.0001", "1e-05", "1e-06"]; tn = dict(zip(TG, "3456"))
RT = "results/L_redteam_gemini.json (code/redteam_checks_gemini.py)"
G = json.load(open(f"{R}/L_redteam_gemini.json"))

# ---------------- (1a) reference-value distance factor, block and per round
A = G["item1a_distance_factor"]
want_pr = {"0.001": 0.98, "0.0001": 1.31, "1e-05": 1.28, "1e-06": 1.08}
for t in TG:
    a = A[t]; n = tn[t]
    assert abs(a["factor_perround"] - want_pr[t]) < 0.006, (t, a["factor_perround"])
    assert abs(a["factor_block"] / a["factor_perround"] - 25 / 21) < 3e-3
    put(f"v6.dist.pr.t{n}", f"{a['factor_perround']:.2f}", RT, "per-round projection of the d=3,5 references to d=7, divided by the d=7 reference")
    put(f"v6.dist.blk.t{n}", f"{a['factor_block']:.2f}", RT, "the same on block rates (registered arithmetic)")
    lf = lambda x: f"{x:.3g}" if x < 1000 else comma(int(float(f"{x:.3g}")))
    put(f"v6.lam35.t{n}", lf(a["lam_pr_35"]), RT, "per-round Lambda, d=3 to 5, reference values (3 significant figures)")
    put(f"v6.lam57.t{n}", lf(a["lam_pr_57"]), RT, "per-round Lambda, d=5 to 7, reference values")
    put(f"v6.lamblk57.t{n}", lf(a["lam_blk_57"]), RT, "block-rate ratio p_L(5)/p_L(7), reference values")
put("v6.dist.pr.min", f"{min(A[t]['factor_perround'] for t in TG):.2f}", RT, ""); put("v6.dist.pr.max", f"{max(A[t]['factor_perround'] for t in TG):.2f}", RT, "")
SMR = G["item2b_smooth_ref"]
for t in ("0.0001", "1e-05"): put(f"v6.dist.sm.t{tn[t]}", f"{A[t]['factor_perround']*A[t]['p7']/SMR[t]:.2f}", RT, "per-round distance factor against the smooth-curve d=7 reference (post hoc)")
put("v6.alg", "25/21", "algebra", "block-rate overshoot of a per-round geometric step from d=3,5 to 7"); put("v6.algval", f"{25/21:.2f}", "algebra", "")
put("v6.alg9", "49/45", "algebra", "the same from d=5,7 to 9")

# ---------------- (1b) R5, block (registered) and per round, both windows
for key, w in (("w1", "item1b_R5_window1"), ("w2", "item1b_R5_bw")):
    B = G[w]; assert B["sanity_mismatch_vs_frozen_csv"] == []
    for t in TG:
        b = B[t]; n = tn[t]; assert b["reps"] == 300
        put(f"v6.r5.{key}.hits.blk.t{n}", b["hits_block"], RT, "R5 hits of 300, block-rate projection (registered)")
        put(f"v6.r5.{key}.hits.pr.t{n}", b["hits_perround"], RT, "R5 hits of 300, per-round projection (corrective)")
        put(f"v6.r5.{key}.cov.blk.t{n}", f"{100*b['cov_block']:.1f}", RT, "percent"); put(f"v6.r5.{key}.cov.pr.t{n}", f"{100*b['cov_perround']:.1f}", RT, "percent")
        f = (lambda x: f"{x:.2f}" if x < 10 else f"{x:.1f}")
        put(f"v6.r5.{key}.rat.blk.t{n}", f(b["med_ratio_block"]), RT, "median R5 / d=7 reference, block"); put(f"v6.r5.{key}.rat.pr.t{n}", f(b["med_ratio_perround"]), RT, "per round")
        put(f"v6.r5.{key}.fit.pr.t{n}", f(b["med_ratio_perround"] / A[t]["factor_perround"]), RT, "per-round R5 median divided by the per-round reference projection (fit-input factor)")
    put(f"v6.r5.{key}.shift", f"{B['0.0001']['med_shift']:.2f}", RT, "median ratio of per-round to block R5 prediction")
W1, W2 = G["item1b_R5_window1"], G["item1b_R5_bw"]
assert [W1[t]["hits_perround"] for t in TG] == [44, 8, 8, 8] and [W2[t]["hits_perround"] for t in TG] == [281, 266, 271, 280]
far2 = [W2[t]["cov_perround"] for t in TG[1:]]
put("v6.r5.w2.far.below90.pr", sum(c < 0.9 for c in far2), RT, "per-round R5 cells at p<=1e-4 below 90%, second window"); put("v6.r5.w2.far.below90.blk", sum(W2[t]["cov_block"] < 0.9 for t in TG[1:]), RT, "")
put("v6.r5.w2.cov.pr.min", f"{100*min(W2[t]['cov_perround'] for t in TG):.0f}", RT, ""); put("v6.r5.w2.cov.pr.max", f"{100*max(W2[t]['cov_perround'] for t in TG):.0f}", RT, "")
put("v6.r5.w2.rat.pr.min", f"{min(W2[t]['med_ratio_perround'] for t in TG):.2f}", RT, ""); put("v6.r5.w2.rat.pr.max", f"{max(W2[t]['med_ratio_perround'] for t in TG):.2f}", RT, "")
put("v6.r5.w1.cov.pr.farmax", f"{100*max(W1[t]['cov_perround'] for t in TG[1:]):.1f}", RT, "")
sh = [math.log(W1[t]["med_ratio_perround"] / A[t]["factor_perround"]) / math.log(W1[t]["med_ratio_perround"]) for t in TG]
put("v6.r5.w1.share.t4", f"{100*sh[1]:.0f}", RT, "share of ln(total per-round R5 error) from the fitted inputs, percent"); put("v6.r5.w1.share.t6", f"{100*sh[3]:.0f}", RT, "")
with open("tab_v6_r5.tex", "w") as f:
    f.write("\\begin{tabular}{lrrrrr}\\hline\\hline\n & \\multicolumn{2}{c}{ref.\\ step / ref.} & \\multicolumn{3}{c}{R5 per round}\\\\\n$p$ & per round & block & med./ref. & fit-input & hits\\\\\\hline\n")
    for t in TG:
        n = tn[t]; f.write(f"$10^{{-{n}}}$ & \\V{{v6.dist.pr.t{n}}} & \\V{{v6.dist.blk.t{n}}} & \\V{{v6.r5.w1.rat.pr.t{n}}} & \\V{{v6.r5.w1.fit.pr.t{n}}} & \\V{{v6.r5.w1.hits.pr.t{n}}}\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")

# ---------------- (1c) d* per round
D = G["item1c_dstar"]; SRC = [("truth", "reference"), ("R1", "R1"), ("R2", "R2"), ("R3", "R3"), ("R4p", "R4$'$"), ("R2lin", "R2-lin"), ("R4", "R4 (diagnostic)"), ("R5", "R5")]
assert {k: D[k]["dstar_perround"] for k, _ in SRC} == dict(truth=21, R1=23, R2=21, R3=21, R4p=21, R2lin=21, R4=21, R5=23)
assert {k: D[k]["dstar_block"] for k, _ in SRC} == dict(truth=23, R1=25, R2=21, R3=23, R4p=23, R2lin=23, R4=23, R5=27) and all(D[k]["dstar_block_matches_frozen"] for k, _ in SRC)
for k, _ in SRC:
    x = D[k]; d = x["dstar_perround"]; assert x["dstar_perround_qubits"] == 2 * d * d and x["dstar_perround_qubit_rounds"] == 2 * d ** 3
    put(f"v6.ds.{k}.lam", f"{x['lam_perround']:.1f}", RT, "per-round Lambda from d=5 to 7 at p=1e-3 (R5: from R1 at d=3 and 5)"); put(f"v6.ds.{k}.d", d, RT, "smallest odd d with block rate <= 1e-12 per d-round experiment, per-round continuation")
    put(f"v6.ds.{k}.q", comma(2 * d * d), RT, "2 d^2"); put(f"v6.ds.{k}.st", comma(2 * d ** 3), RT, "2 d^3"); put(f"v6.ds.{k}.dblk", x["dstar_block"], RT, "d* by the block-rate rule (registered arithmetic)")
put("v6.ds.inflation", f"{D['R1']['dstar_perround_qubit_rounds']/D['truth']['dstar_perround_qubit_rounds']:.2f}", RT, "qubit-rounds of R1 and R5 over the reference, per round")
with open("tab_v6_dstar.tex", "w") as f:
    f.write("\\begin{tabular}{lrrrrr}\\hline\\hline\nsource & $\\pL(7)$/ref. & $\\Lambda$ & $d^*$ & qubits & qubit-rounds\\\\\\hline\n")
    for k, name in SRC: f.write(f"{name} & \\V{{am.ds.{k}.rel}} & \\V{{v6.ds.{k}.lam}} & \\V{{v6.ds.{k}.d}} & \\V{{v6.ds.{k}.q}} & \\V{{v6.ds.{k}.st}}\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")

# ---------------- (2) the d=7 reference: combination, shape, R2 re-scored
C = G["item2a_combination"]
put("v6.eq.maxdiff", f"{100*max(abs(C[t]['eq_over_ivw']) for t in TG):.1f}", RT, "largest |equal-weight / inverse-variance - 1| over the four d=7 cells, percent")
put("v6.run.zmax", f"{max(abs(C[t]['z_run_diff']) for t in TG):.1f}", RT, "largest |z| between the two runs of a d=7 cell")
for t in TG:
    n = tn[t]; s = float(t) ** 4 * 1e6
    for lab in ("ivw", "run1", "run2"): put(f"v6.d7c.{lab}.t{n}", f"{C[t][lab]/s:.2f}", RT, "p_L / p^4 in units of 1e6, d=7")
    put(f"v6.d7c.w1.t{n}", f"{C[t]['weight1']:.2f}", RT, "inverse-variance weight of run 1"); put(f"v6.d7c.z.t{n}", f"{C[t]['z_run_diff']:+.2f}", RT, "z of run 1 minus run 2")
S = G["item2b_shape"]; assert abs(S["model_A_low"]["chi2"] - 6.42) < 0.01 and S["model_A_low"]["dof"] == 5 and S["model_A_B"]["dof"] == 6
put("v6.smooth.chi2", f"{S['model_A_low']['chi2']:.1f}", RT, "A exp(-lambda) on the six runs at p<=1e-4 (post hoc fallback)"); put("v6.smooth.dof", 5, RT, ""); put("v6.smooth.p", f"{S['model_A_low']['p_value']:.2f}", RT, "")
put("v6.smooth.A", f"{S['model_A_low']['A']:.2f}", RT, "fitted A in units of 1e6"); put("v6.smooth.zmax", f"{max(abs(z) for z in S['model_A_low']['z']):.1f}", RT, "")
put("v6.smooth.reg.chi2", f"{S['model_A_B']['chi2']:.1f}", RT, "pre-registered model A exp(-lambda)(1+Bp), all 8 runs: failed"); put("v6.smooth.reg.dof", 6, RT, "")
put("v6.dip.se.min", "7", "numbers.json truthse.d7.*", ""); put("v6.dip.se.max", "8", "numbers.json truthse.d7.*", "")
VAR = [("ivw", "registered"), ("equal", "equal weights"), ("run1", "run 1"), ("run2", "run 2"), ("smoothed_model", "smooth curve")]
frozen = {}
for fn, key in (("L_coverage.csv", "w1"), ("L_bw_coverage.csv", "w2")):
    for r in csv.DictReader(open(f"{R}/{fn}")): frozen[(key, r["recipe"], int(r["d"]), float(r["target"]))] = int(r["hits"])
for key, w in (("w1", "item2c_R2_d7_window1"), ("w2", "item2c_R2_d7_bw")):
    H = G[w]; assert H["sanity_mismatch"] == [] and H["reps"] == 300
    for t in TG:
        n = tn[t]; assert H["hits"]["ivw"][t] == frozen[(key, "R2", 7, float(t))]
        hs = [H["hits"][v][t] for v, _ in VAR]
        for v, _ in VAR: put(f"v6.r2d7.{key}.{v.split('_')[0]}.t{n}", f"{100*H['hits'][v][t]/300:.0f}", RT, "R2 d=7 hit fraction (percent of 300) against this reference")
        put(f"v6.r2d7.{key}.min.t{n}", f"{100*min(hs)/300:.0f}", RT, "smallest over the five references"); put(f"v6.r2d7.{key}.max.t{n}", f"{100*max(hs)/300:.0f}", RT, "largest over the five references")
    rr = [H["med_ratio"][v][t] for v, _ in VAR for t in TG[1:3]]
    put(f"v6.r2d7.{key}.rat.min", f"{min(rr):.2f}", RT, "R2 d=7 median prediction / reference at 1e-4 and 1e-5, smallest over references"); put(f"v6.r2d7.{key}.rat.max", f"{max(rr):.2f}", RT, "largest")
assert V["v6.r2d7.w1.max.t4"][0] == "92" and V["v6.r2d7.w1.min.t4"][0] == "0"
put("v6.r2d7.eq.maxmove", max(abs(G[w]["hits"]["equal"][t] - G[w]["hits"]["ivw"][t]) for w in ("item2c_R2_d7_window1", "item2c_R2_d7_bw") for t in TG), RT, "largest change in hits of 300 from equal weights")
with open("tab_v6_r2d7.tex", "w") as f:
    f.write("\\begin{tabular}{lrrrrrrrr}\\hline\\hline\n & \\multicolumn{4}{c}{first window} & \\multicolumn{4}{c}{second window}\\\\\nreference & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for v, name in VAR:
        f.write(name + "".join(f" & \\V{{v6.r2d7.{key}.{v.split('_')[0]}.t{tn[t]}}}" for key in ("w1", "w2") for t in TG) + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
with open("tab_v6_d7runs.tex", "w") as f:
    f.write("\\begin{tabular}{lrrrr}\\hline\\hline\n$\\pL/p^4$ ($10^6$), $d=7$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for lab, name in (("run1", "run 1"), ("run2", "run 2"), ("ivw", "combined (frozen)")): f.write(name + "".join(f" & \\V{{v6.d7c.{lab}.t{tn[t]}}}" for t in TG) + "\\\\\n")
    f.write("weight of run 1" + "".join(f" & \\V{{v6.d7c.w1.t{tn[t]}}}" for t in TG) + "\\\\\n$z$ (run 1 $-$ run 2)" + "".join(f" & $\\V{{v6.d7c.z.t{tn[t]}}}$" for t in TG) + "\\\\\n\\hline\\hline\n\\end{tabular}\n")

# ---------------- (3) pre-registered, not yet run: d=7 refinement, sliding sweep
P3 = G["item3_plan"]; rate = P3["_rate_draws_per_s_per_core"]["median"]
tot = sum(P3[t]["Mtot_plain_1p5pct"] for t in TG); tots = sum(P3[t]["Mtot_plain_x1p5_safety"] for t in TG)
h1, h2 = tot / rate / 3600, tots / rate / 3600; assert 95 < h1 < 105 and 143 < h2 < 157, (h1, h2)
put("v6.ref7.cpuh", "100", RT, f"draws for 1.5% plain se, summed over the four d=7 cells, at the median measured rate: {h1:.1f} CPU-h"); put("v6.ref7.cpuh.safe", "150", RT, f"with the x1.5 safety factor: {h2:.1f} CPU-h")
put("v6.ref7.draws", sci(tots, 2), RT, "draws with the safety factor"); put("v6.ref7.bar", "1.5", "PLAN.md amendment item 3", "target relative se, percent")
SW = G["itema_sweep_cost"]
put("v6.sweep.cpuh.min", "6", RT, f"windows at the measured rate: {SW['d3']['total_cpu_h_measured']+SW['d5']['total_cpu_h_measured']:.1f} CPU-h"); put("v6.sweep.cpuh.max", "35", "PLAN.md addendum (a)", "windows at a third of the rate plus references")

# ---------------- crossings, tail, surrogate bracket, S16
X = G["item5_crossings_1e-3"]
for a, b in ((3, 5), (3, 7), (5, 7)):
    put(f"v6.cross.{a}{b}", f"{X[f'mean_d{a}_d{b}'][0]:.2f}", RT, "crossing of the window means, units of 1e-3"); put(f"v6.cross.med.{a}{b}", f"{X[f'median_d{a}_d{b}'][0]:.2f}", RT, "crossing of the window medians")
T = G["itemc_tail120"]
put("v6.tail.d9", f"10^{{{math.ceil(T['d9_p0.0001']['chernoff_log10_bound_P_W_ge_121'])}}}", RT, "Chernoff bound on P(W>=121), d=9, p=1e-4"); put("v6.tail.d7", f"10^{{{math.ceil(T['d7_p0.001']['chernoff_log10_bound_P_W_ge_121'])}}}", RT, "d=7, p=1e-3 (largest mean number of faults)")
BR = G["itemb_surrogate_bracket"]; ins = tot_r = isr = ist = 0
with open("tab_v6_bracket.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrrr}\\hline\\hline\n$d$ & $p$ & 4th-order polynomial & certified $[L,U]$, full value & orders 5--7 / poly. & $(U-L)$/poly. & intervals containing $[L,U]$\\\\\\hline\n")
    for d in (3, 5):
        for t in ("0.001", "0.0001"):
            b = BR[f"d{d}_p{t}"]; k = f"v6.br.d{d}.t{tn[t]}"; w = b["rel_width"]
            nd = 7 if w < 1e-4 else 4; kn = 100 * b["known_higher_order_over_poly4"]
            put(k + ".poly", sci(b["poly4"], nd), RT, "surrogate polynomial through order 4"); put(k + ".L", sci(b["L_exact"], nd), RT, "certified lower end, full surrogate value"); put(k + ".U", sci(b["U_exact"], nd), RT, "certified upper end")
            put(k + ".known", ("-" + sci(-kn, 1)) if abs(kn) < 1e-3 else f"{kn:.2g}", RT, "orders 5-7 of shells <=4 relative to the polynomial, percent")
            put(k + ".width", sci(100 * w, 2) if w < 1e-4 else (f"{100*w:.2g}" if w < 10 else f"{round(100*w, -2):.0f}"), RT, "bracket width (U-L) over polynomial, always percent")
            put(k + ".in", b["bracket_inside"], RT, "certified replicate intervals (of 24) containing the whole bracket"); put(k + ".n", b["n_reps"], RT, "")
            assert b["inside_poly4"] == b["n_reps"] == 24 and b["bracket_intersects"] == 24 and b["is_runs"] == 20
            ins += b["inside_poly4"]; tot_r += b["n_reps"]; isr += b["is_inside_poly4"]; ist += b["is_touch_bracket"]
            f.write(f"{d} & $10^{{-{tn[t]}}}$ & $\\V{{{k}.poly}}$ & $\\V{{{k}.L}}$ to $\\V{{{k}.U}}$ & $\\V{{{k}.known}}$\\% & $\\V{{{k}.width}}$\\% & \\V{{{k}.in}}/\\V{{{k}.n}}\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
assert (ins, tot_r, isr, ist) == (96, 96, 78, 78)
put("v6.br.is.touch", ist, RT, "IS runs (of 80) whose interval touches the certified bracket"); put("v6.br.is.n", 80, RT, "")
put("v6.br.d3.width.max", f"{100*BR['d3_p0.001']['rel_width']:.1f}", RT, "bracket width at d=3, p=1e-3, percent"); put("v6.br.d5.t4.pct", f"{100*BR['d5_p0.0001']['rel_width']:.0f}", RT, "percent")
put("v6.br.d5.t3.lower", f"{-100*BR['d5_p0.001']['known_higher_order_over_poly4']:.0f}", RT, "known orders 5-7 lower the polynomial by this percent at d=5, p=1e-3")
GR = G["itemc_grids"]
gl = [k for k in GR if isinstance(GR[k], list) and len(GR[k]) == 18 and "admissible" not in k]; assert len(gl) == 2, gl
for k in gl:
    lab = "upper" if max(GR[k]) < 1 else "lower"
    put(f"v6.grid.{lab}", ", ".join(f"{x:.3g}" for x in GR[k]), RT, f"the 18 fixed multipliers g_j (lambda_j = g_j / b_w), {lab} limit")
put("v6.grid.lower.full", sum(x == 1.0 for x in GR["admissible_m_over_b_for_Kplus"]), RT, "number of lower-limit multipliers admissible for every trial mean")
S16 = G["item4_S16"]
put("v6.s16.ratio11", f"{S16['strata']['11']['ucb_over_all_zero']-1:.0e}".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}", RT, "relative excess of the w=11 upper limit over the zero-data limit")
put("v6.s16.f11", f"{S16['strata']['11']['ucb_actual']:.4f}", RT, "upper limit on f_11"); put("v6.s16.f12", f"{S16['strata']['12']['ucb_actual']:.4f}", RT, "upper limit on f_12")
put("v6.s16.ymax", sci(S16["strata"]["12"]["Ymax"], 2), RT, "largest score of a failing draw in strata 11 and 12"); assert S16["strata"]["12"]["Ymax"] > S16["strata"]["11"]["Ymax"]

# ---------------- coverage with the reference range, both windows
AM = json.load(open(f"{R}/L_amend_gemini.json"))["reference_uncertainty"]; BW = json.load(open(f"{R}/L_r4p_bw.json")); assert BW["_sanity_mismatch"] == 0
RU = "results/L_amend_gemini.json (first window), results/L_r4p_bw.json (second window)"
cell = {}
for rc in ("R1", "R2", "R3", "R4", "R4p", "R5", "R2lin"):
    for d in ((7,) if rc == "R5" else (3, 5, 7)):
        for t in TG:
            k = f"{rc}|{d}|{t}|primary"; a = AM[k]; b = BW[k]
            if rc in ("R1", "R2", "R3", "R4", "R5"):
                assert round(a["hit_at_point"] * 300) == frozen[("w1", rc, d, float(t))] and b["hits"] == frozen[("w2", rc, d, float(t))], k
            cell[("w1", rc, d, t)] = (a["hit_at_point"], a["hit_min"], a["hit_max"]); cell[("w2", rc, d, t)] = (b["cov"], b["hit_min"], b["hit_max"])
for (w, rc, d, t), (c, lo, hi) in cell.items():
    if rc in ("R4", "R2lin"): continue
    put(f"v6.cv.{w}.{rc}.d{d}.t{tn[t]}", f"{100*c:.1f}", RU, "coverage at the reference point estimate, percent")
    put(f"v6.rg.{w}.{rc}.d{d}.t{tn[t]}", f"{100*lo:.0f}--{100*hi:.0f}", RU, "smallest and largest hit fraction as one shared reference moves over estimate +/- 1.96 se, percent")
def count(w, rcs, pred, ds=(3, 5, 7), ts=TG): return sum(pred(*cell[(w, rc, d, t)]) for rc in rcs for d in ds for t in ts if (w, rc, d, t) in cell)
below = lambda c, lo, hi: hi < 0.90          # below 90% for every reference value in the range
strad = lambda c, lo, hi: lo < 0.90 <= hi    # the 90% line is crossed inside the range
above = lambda c, lo, hi: lo >= 0.90
for w in ("w1", "w2"):
    fx = ("R2", "R3")
    put(f"v6.fx.{w}.pt", count(w, fx, lambda c, lo, hi: c < 0.90), RU, "R2+R3 cells (of 24) below 90% at the point estimate")
    put(f"v6.fx.{w}.robust", count(w, fx, below), RU, "R2+R3 cells (of 24) below 90% at every reference value in the range")
    put(f"v6.fx.{w}.strad", count(w, fx, strad), RU, "R2+R3 cells whose range crosses 90%"); put(f"v6.fx.{w}.above", count(w, fx, above), RU, "R2+R3 cells at or above 90% over the whole range")
    put(f"v6.fx.{w}.robust.lowd", count(w, fx, below, ds=(3, 5)), RU, "of 16 at d=3,5"); put(f"v6.fx.{w}.robust.d7", count(w, fx, below, ds=(7,)), RU, "of 8 at d=7")
    put(f"v6.fx.{w}.strad.d7", count(w, fx, strad, ds=(7,)), RU, "of 8 at d=7"); put(f"v6.fx.{w}.strad.lowd", count(w, fx, strad, ds=(3, 5)), RU, "of 16 at d=3,5")
    put(f"v6.R1.{w}.robust", count(w, ("R1",), below), RU, "R1 cells (of 12) below 90% over the whole range"); put(f"v6.R1.{w}.above", count(w, ("R1",), above), RU, "R1 cells (of 12) at or above 90% over the whole range")
    put(f"v6.R1.{w}.strad", count(w, ("R1",), strad), RU, "")
    put(f"v6.R4p.{w}.lowd.strad", count(w, ("R4p",), strad, ds=(3, 5)), RU, "R4' cells at d=3,5 whose range crosses 90%"); put(f"v6.R4p.{w}.d7.strad", count(w, ("R4p",), strad, ds=(7,)), RU, "R4' cells at d=7 whose range crosses 90%")
    lows = [cell[(w, "R4p", d, t)] for d in (3, 5) for t in TG]
    put(f"v6.R4p.{w}.lowd.rgmin", f"{100*min(x[1] for x in lows):.0f}", RU, "smallest hit fraction in the range, R4' at d=3,5"); put(f"v6.R4p.{w}.lowd.rgmax", f"{100*max(x[2] for x in lows):.0f}", RU, "largest")
    d7 = [cell[(w, "R4p", 7, t)] for t in TG]
    put(f"v6.R4p.{w}.d7.rgmin", f"{100*min(x[1] for x in d7):.0f}", RU, ""); put(f"v6.R4p.{w}.d7.rgmax", f"{100*max(x[2] for x in d7):.0f}", RU, "")
    r1 = [cell[(w, "R1", d, t)] for d in (3, 5, 7) for t in TG]
    put(f"v6.R1.{w}.rgmin", f"{100*min(x[1] for x in r1):.0f}", RU, ""); put(f"v6.R1.{w}.rgmax", f"{100*max(x[2] for x in r1):.0f}", RU, "")
    r5 = [cell[(w, "R5", 7, t)] for t in TG]
    put(f"v6.R5.{w}.rgmin", f"{100*min(x[1] for x in r5):.0f}", RU, "registered (block-rate) R5"); put(f"v6.R5.{w}.rgmax", f"{100*max(x[2] for x in r5):.0f}", RU, "")
ROWS = [(rc, d) for rc in ("R2", "R3", "R1", "R4p") for d in (3, 5, 7)] + [("R5", 7)]; NAME = dict(R1="R1", R2="R2", R3="R3", R4p="R4$'$", R5="R5")
with open("tab_v6_cov.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrrrrrr}\\hline\\hline\n & & \\multicolumn{4}{c}{first window ($4$--$15\\times10^{-3}$)} & \\multicolumn{4}{c}{second window ($1.5$--$4\\times10^{-3}$)}\\\\\nrecipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for rc, d in ROWS:
        f.write(f"{NAME[rc]} & {d}" + "".join(f" & \\V{{v6.cv.{w}.{rc}.d{d}.t{tn[t]}}} (\\V{{v6.rg.{w}.{rc}.d{d}.t{tn[t]}}})" for w in ("w1", "w2") for t in TG) + "\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
def rk(w, rc, d, n):
    if rc == "R4p": return f"am.ratio.R4p.d{d}.t{n}" if w == "w1" else f"r4b.ratio.d{d}.t{n}"
    return f"ratio.{rc}.d{d}.t{n}" if w == "w1" else f"bw.ratio.{rc}.d{d}.t{n}"
with open("tab_v6_ratio.tex", "w") as f:
    f.write("\\begin{tabular}{llrrrrrrrr}\\hline\\hline\n & & \\multicolumn{4}{c}{first window} & \\multicolumn{4}{c}{second window}\\\\\nrecipe & $d$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$ & $10^{-3}$ & $10^{-4}$ & $10^{-5}$ & $10^{-6}$\\\\\\hline\n")
    for rc, d in ROWS:
        f.write(f"{NAME[rc]}{' (block)' if rc == 'R5' else ''} & {d}" + "".join(f" & \\V{{{rk(w, rc, d, tn[t])}}}" for w in ("w1", "w2") for t in TG) + "\\\\\n")
    f.write("R5 (per round) & 7" + "".join(f" & \\V{{v6.r5.{w}.rat.pr.t{tn[t]}}}" for w in ("w1", "w2") for t in TG) + "\\\\\n\\hline\\hline\n\\end{tabular}\n")

with open("numbers_v6.tex", "w") as f:
    f.write("% GENERATED by make_numbers_v6.py. Do not edit. Input after numbers_v5.tex.\n")
    for k, (v, s, h) in V.items(): f.write(f"\\expandafter\\def\\csname LV@{k}\\endcsname{{{v}}}\n")
json.dump({k: dict(value=v, source=s, how=h) for k, (v, s, h) in V.items()}, open("numbers_v6.json", "w"), indent=1)
for k, (v, s, h) in V.items():
    if not k.startswith(("v6.cv.", "v6.rg.")): print(k, v)
