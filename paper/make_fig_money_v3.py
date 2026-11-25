"""Fig. 1 of draft v3 (plotting only, no inference): the window, the reference estimates and the reported intervals of R1, R2, R4 and R4'.
From stored files only. Run from this folder:  python3 make_fig_money_v3.py
Top row: p_L against p. Bottom row: the same divided by p^((d+1)/2), each panel on its own narrow axis so the prefactor drift is visible.
Ribbons are the medians over the 300 replicates of the lower and upper 95% limits each recipe reports (primary analysis)."""
import sys, json, math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../code")
from coverage import fit_all, load_reps, wls, PS, X0
import os; os.chdir("../code"); from analyze_main import load_truth; T = load_truth(); os.chdir("../paper")
R = "../results"; DS = (3, 5, 7); reps = load_reps(f"{R}/fit_d357.json"); assert len(reps) == 300
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "font.family": "serif",
                     "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
COL = {"R1": "#0072B2", "R2": "#D55E00", "R4": "#CC79A7", "R4p": "#009E73"}; LAB = {"R1": "R1", "R2": "R2", "R4": "R4 (registered)", "R4p": "R4$'$ (corrective)"}
def r4p(data, tp):   # identical to R4 in coverage.fit_all except the exponent (as code/amend_gemini.py)
    out = {}
    for d in data:
        k, n = data[d]; P = k / n; h = (d + 1) / 2; X = np.c_[np.ones(5), PS / 1e-2, (PS / 1e-2) ** 2]
        b, C, c2, dof = wls(X, np.log(P) - h * X0, k / (1 - P)); g = np.array([1, tp / 1e-2, (tp / 1e-2) ** 2]); out[("R4p", d)] = (g @ b + h * np.log(tp), g @ C @ g)
    return out
pg = np.unique(np.concatenate([np.logspace(-6, math.log10(1.5e-2), 90), [1e-3, 1e-4, 1e-5, 1e-6]])); band = {}
for tp in pg:
    acc = {}
    for rep in reps:
        o = fit_all({d: rep[d] for d in DS}, tp); o.update(r4p({d: rep[d] for d in DS}, tp))
        for (rc, d), (mu, var) in o.items():
            if rc in COL: acc.setdefault((rc, d), []).append((mu - 1.96 * math.sqrt(var), mu, mu + 1.96 * math.sqrt(var)))
    for k, v in acc.items(): band.setdefault(k, []).append(np.median(np.array(v), 0))
band = {k: np.exp(np.array(v)) for k, v in band.items()}
# sanity: medians at the targets reproduce the committed numbers
NA = json.load(open("numbers_amend.json")); NN = json.load(open("numbers.json"))
for d in DS:
    for t, tn in ((1e-3, 3), (1e-4, 4), (1e-5, 5), (1e-6, 6)):
        i = int(np.argmin(abs(pg - t))); ref = T[(d, t)]["truth"]
        assert f"{band[('R4p', d)][i, 1] / ref:.2f}" == NA[f"am.ratio.R4p.d{d}.t{tn}"]["value"], (d, t)
        for rc in ("R1", "R4"):
            x = band[(rc, d)][i, 1] / ref; want = float(NN[f"ratio.{rc}.d{d}.t{tn}"]["value"]); assert abs(x / want - 1) < 0.06, (rc, d, t, x, want)
W = {d: np.median(np.array([r[d][0] / r[d][1] for r in reps]), 0) for d in DS}
ld = np.log(W[5]) - np.log(W[7]); cr = [math.exp(math.log(PS[i]) + ld[i] / (ld[i] - ld[i + 1]) * (math.log(PS[i + 1]) - math.log(PS[i]))) for i in range(4) if ld[i] * ld[i + 1] < 0]
assert len(cr) == 1 and f"{cr[0]*1e3:.1f}" == NA["am.cross"]["value"]
EX = {}
for d in (3, 5):
    c = json.load(open(f"../exact/d{d}_manifest.json"))["coefficients_p_L=sum c_k p^k"]; assert all(v["status"] == "EXACT" for v in c.values())
    EX[d] = {int(k): v["coef"] for k, v in c.items()}
def draw(row, fname, figsize):
    """row 0: p_L against p (SM figure); row 1: p_L / p^((d+1)/2) (main-text figure). Same data, same code path as the original two-row figure."""
    fig, axs = plt.subplots(1, 3, figsize=figsize, gridspec_kw=dict(wspace=0.28))
    nrm0 = (lambda h: (lambda p: 1.0)) if row == 0 else (lambda h: (lambda p: p ** h))
    for j, d in enumerate(DS):
        h = (d + 1) / 2; nrm = nrm0(h); ax = axs[j]; tp = sorted((p, x["truth"], x["se"]) for (dd, p), x in T.items() if dd == d)
        P = np.array([a[0] for a in tp]); Y = np.array([a[1] for a in tp]); E = 1.96 * np.array([a[2] for a in tp])
        for rc in (("R1", "R4", "R2", "R4p") if row == 0 else ("R1", "R2", "R4p")):   # v3 (arXiv): R4 (wrong exponent) is drawn only in the SM figure
            b = band[(rc, d)]; ax.fill_between(pg, b[:, 0] / nrm(pg), b[:, 2] / nrm(pg), color=COL[rc], alpha=.28, lw=0)
            ax.plot(pg, b[:, 1] / nrm(pg), color=COL[rc], lw=.9, label=LAB[rc] if j == 0 else None)
        if row == 1 and d in EX:
            ax.plot(pg, sum(c * pg ** (k - h) for k, c in EX[d].items()), ":", color="k", lw=1.1, zorder=6, label="surrogate, through $p^4$" if d == 3 else None)
        ax.errorbar(P, Y / nrm(P), yerr=E / nrm(P), fmt="o", ms=3, color="k", lw=.8, zorder=5, label="reference estimates" if j == 0 else None)
        ax.plot(PS, W[d] / nrm(PS), "s", mfc="none", ms=5, color="k", mew=.8, zorder=4, label="window medians" if j == 0 else None)
        ax.axvspan(4e-3, 1.5e-2, color="gray", alpha=.10, lw=0); ax.axvline(cr[0], color="gray", ls="--", lw=.7)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(7e-7, 2e-2); ax.set_xlabel("physical error rate $p$")
        if row == 0: ax.set_ylim(1e-21, 1); ax.set_title(f"$d={d}$", fontsize=8)
        else:
            plateau = T[(d, 1e-6)]["truth"] / 1e-6 ** h; ax.set_ylim(plateau / 2.5, plateau * 3.2); ax.set_title(f"$d={d}$", fontsize=8)
    axs[0].set_ylabel("$p_L$" if row == 0 else "$p_L/p^{(d+1)/2}$")
    if row == 0: axs[0].legend(loc="lower right", frameon=False, handlelength=1.2, fontsize=6, borderaxespad=0.1, labelspacing=0.3)
    else:
        hs, ls = axs[0].get_legend_handles_labels(); fig.legend(hs, ls, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False, handlelength=1.4, fontsize=6.5, columnspacing=1.0)
    fig.savefig(fname + ".pdf", bbox_inches="tight"); fig.savefig(fname + ".png", bbox_inches="tight", dpi=200); plt.close(fig)
# v9 (arXiv v2): the lower row stays in the main text as figs/fig0_money, the upper row moves to the supplement as figs/fig0_money_top
draw(1, "figs/fig0_money", (7.1, 2.5)); draw(0, "figs/fig0_money_top", (7.1, 2.7)); print("ok", cr)
