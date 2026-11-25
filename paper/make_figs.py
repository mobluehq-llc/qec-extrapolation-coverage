"""Figures for paper L, from the committed result files only. Run from this folder:
python3 make_figs.py"""
import sys, json, csv, glob, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "../code")
from coverage import fit_all, load_reps
R = "../results"; EX = "../exact"; O = "figs"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 300})
COL = {"R1": "#0072B2", "R2": "#D55E00", "R3": "#009E73", "R4": "#CC79A7", "R5": "#000000"}
MK = {3: "o", 5: "s", 7: "^"}
RC = ["R1", "R2", "R3", "R4", "R5"]; DS = [3, 5, 7]; TG = [1e-3, 1e-4, 1e-5, 1e-6]
def save(fig, n):
    fig.savefig(f"{O}/{n}.pdf", bbox_inches="tight"); fig.savefig(f"{O}/{n}.png", bbox_inches="tight", dpi=200); plt.close(fig)
T = {}
for f in glob.glob(f"{R}/truth_d*_p*.json"):
    x = json.load(open(f)); T[(x["d"], x["p"])] = x
for f in glob.glob(f"{R}/truthextra_d*_p*.json"):
    y = json.load(open(f)); x = T[(y["d"], y["p"])]; w1, w2 = 1 / x["se"] ** 2, 1 / y["se"] ** 2
    x["truth"], x["se"] = (w1 * x["truth"] + w2 * y["truth"]) / (w1 + w2), (w1 + w2) ** -0.5
truth = {(d, p): T[(d, p)]["truth"] for d in DS for p in TG}

# ---- Fig. 1: coverage heat-map, recipe x d (rows) by decades (columns); primary and sensitivity
rows = [(rc, d) for rc in RC for d in DS if not (rc == "R5" and d != 7)]
fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.9), sharey=True)
for ax, fn, ttl in zip(axs, ("L_coverage.csv", "L_coverage_scaled.csv"), ("(a) primary: absolute $\\sigma$", "(b) sensitivity: $\\chi^2$/dof-scaled")):
    cs = {(r["recipe"], int(r["d"]), float(r["target"])): float(r["cov"]) for r in csv.DictReader(open(f"{R}/{fn}"))}
    A = np.array([[100 * cs[(rc, d, p)] for p in TG] for rc, d in rows])
    im = ax.imshow(A, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    for i in range(A.shape[0]):
        for j in range(4):
            ax.text(j, i, f"{A[i, j]:.0f}" if A[i, j] in (0, 100) else f"{A[i, j]:.1f}", ha="center", va="center", fontsize=6.5, color="white" if A[i, j] > 55 else "black")
    for y in (2.5, 5.5, 8.5, 11.5): ax.axhline(y, color="white", lw=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(["$10^{-3}$\n(1 decade)", "$10^{-4}$\n(2)", "$10^{-5}$\n(3)", "$10^{-6}$\n(4)"])
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{rc}, $d={d}$" for rc, d in rows])
    ax.set_xlabel("target physical error rate $p$"); ax.set_title(ttl, fontsize=8)
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
cb = fig.colorbar(im, ax=axs, fraction=0.03, pad=0.02); cb.set_label("coverage of the reference estimate by the 95% interval (%); nominal 95")
save(fig, "fig1_coverage")

# ---- Fig. 2: predicted against true logical error rate
reps = load_reps(f"{R}/fit_d357.json"); P = {}
for data in reps:
    for p in TG:
        for (rc, d), (mu, var) in fit_all({d: data[d] for d in DS}, p).items(): P.setdefault((rc, d, p), []).append(mu)
fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.2))
ax = axs[0]
ax.plot([1e-19, 1e-2], [1e-19, 1e-2], color="#888888", lw=0.8, zorder=0)
for (rc, d, p), v in P.items():
    v = np.exp(np.array(v)); lo, md, hi = np.percentile(v, [2.5, 50, 97.5]); t = truth[(d, p)]
    ax.errorbar([t], [md], yerr=[[md - lo], [hi - md]], fmt=MK[d], color=COL[rc], ms=3.5, lw=0.7, mfc="none" if rc in ("R2", "R3") else COL[rc], mew=0.8)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(2e-19, 2e-3); ax.set_ylim(2e-19, 2e-3)
ax.set_xlabel("reference estimate of $p_L$ (real circuit)"); ax.set_ylabel("predicted $p_L$ (median, central 95% of replicates)"); ax.set_title("(a)", fontsize=8, loc="left")
ax.legend(handles=[Line2D([], [], color=COL[r], marker="o", ls="", ms=4, mfc="none" if r in ("R2", "R3") else COL[r], label=r) for r in RC] +
          [Line2D([], [], color="k", marker=MK[d], ls="", ms=4, mfc="none", label=f"$d={d}$") for d in DS], frameon=False, ncol=2, loc="upper left", handletextpad=0.2, columnspacing=0.8)
ax = axs[1]
ax.axhline(1, color="#888888", lw=0.8)
for rc in RC:
    for d in DS:
        if (rc, d, 1e-3) not in P: continue
        r = [float(np.exp(np.median(P[(rc, d, p)]))) / truth[(d, p)] for p in TG]
        ax.plot([1, 2, 3, 4], r, marker=MK[d], color=COL[rc], ms=3.5, lw=0.8, mfc="none" if rc in ("R2", "R3") else COL[rc], mew=0.8)
ax.set_yscale("log"); ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel("decades of extrapolation below $p=10^{-2}$"); ax.set_ylabel("median prediction / reference"); ax.set_title("(b)", fontsize=8, loc="left")
ax.set_ylim(0.5, 150)
fig.tight_layout(); save(fig, "fig2_pred_vs_true")

# ---- Fig. 3: exact low-order shells of the linearised-DEM surrogate, and surrogate against real circuit
M = {d: json.load(open(f"{EX}/d{d}_manifest.json")) for d in DS}; S7 = json.load(open(glob.glob(f"{EX}/d7_w4_snapshot_*.json")[0]))
CD = {3: "#0072B2", 5: "#D55E00", 7: "#009E73"}
fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.9))
ax = axs[0]; FLOOR = 1e-12
for d in DS:
    ws = sorted(int(w) for w in M[d]["shells"]); fr = [M[d]["shells"][str(w)]["n_fail"] / M[d]["shells"][str(w)]["sets_total"] for w in ws]
    nz = [(w, f) for w, f in zip(ws, fr) if f > 0]; z = [w for w, f in zip(ws, fr) if f == 0]
    off = {3: -0.08, 5: 0, 7: 0.08}[d]
    if d == 7:
        e = S7["ESTIMATE_full_shell_NOT_exact"]["n_fail"]; f4 = e["estimate"] / S7["sets_total"]; s4 = e["se"] / S7["sets_total"]
        ax.errorbar([4 + off], [f4], yerr=[1.96 * s4], fmt=MK[d], color=CD[d], ms=4, mfc="none", lw=0.8)
    if nz: ax.plot([w + off for w, f in nz], [f for w, f in nz], marker=MK[d], color=CD[d], ms=4, lw=0.8, label=f"$d={d}$")
    else: ax.plot([], [], marker=MK[d], color=CD[d], ms=4, lw=0.8, label=f"$d={d}$")
    ax.plot([w + off for w in z], [FLOOR] * len(z), marker=MK[d], color=CD[d], ms=4, ls="", clip_on=False)
ax.set_yscale("log"); ax.set_ylim(FLOOR, 1); ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel("number of faults in the set (weight)")
ax.set_ylabel("fraction of fault sets that fail"); ax.set_title("(a)", fontsize=8, loc="left")
ax.text(1.0, 2.5e-12, "exactly zero", fontsize=6.5, color="#555555", va="bottom")
ax.text(3.6, 2e-5, "sampled", fontsize=6.5, color=CD[7], ha="center")
ax.legend(frameon=False, loc="upper left")
ax = axs[1]; ax.axhline(1, color="#888888", lw=0.8)
for d in (3, 5):
    c = {int(k): v["coef"] for k, v in M[d]["coefficients_p_L=sum c_k p^k"].items()}
    r = [truth[(d, p)] / sum(c[k] * p ** k for k in c) for p in TG]; e = [1.96 * T[(d, p)]["se"] / truth[(d, p)] * rr for p, rr in zip(TG, r)]
    ax.errorbar([p * (1.12 if d == 5 else 0.89) for p in TG], r, yerr=e, fmt=MK[d], color=CD[d], ms=4, lw=0.8, label=f"$d={d}$")
ax.set_xscale("log"); ax.set_ylim(0.9, 1.1); ax.set_xlabel("physical error rate $p$"); ax.set_ylabel("real-circuit reference / surrogate polynomial"); ax.set_title("(b)", fontsize=8, loc="left")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout(); save(fig, "fig3_exact_shells")
print("figures written")
