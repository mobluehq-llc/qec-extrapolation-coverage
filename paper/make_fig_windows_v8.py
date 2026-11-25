"""Fig. 2 of v8 (plotting only; replaces the figure made inside make_numbers_v4.py, whose 6.6 x 2.7 in layout was shrunk to column width and left its legend over the data).
Reads results/L_coverage.csv, results/L_bw_coverage.csv, numbers.json (ratio.R1.*), numbers_v4.json (bw.ratio.R1.*). Run from this folder on a host with matplotlib."""
import csv, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R = "../results"; DS = (3, 5, 7); TG = ["0.001", "0.0001", "1e-05", "1e-06"]; XT = [1e-3, 1e-4, 1e-5, 1e-6]
def cov(f): return {(int(r["d"]), r["target"]): 100 * float(r["cov"]) for r in csv.DictReader(open(f"{R}/{f}")) if r["recipe"] == "R1"}
c1, c2 = cov("L_coverage.csv"), cov("L_bw_coverage.csv"); n0 = json.load(open("numbers.json")); n4 = json.load(open("numbers_v4.json"))
tn = dict(zip(TG, "3456")); g = lambda D, k: float(D[k]["value"])
plt.rcParams.update({"font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.0)); MK = {3: "o", 5: "s", 7: "^"}; B, O = "#0072B2", "#D55E00"
for d in DS:
    ax[0].plot(XT, [c1[(d, t)] for t in TG], MK[d] + "--", color=B, mfc="none", ms=3.2, lw=0.8); ax[0].plot(XT, [c2[(d, t)] for t in TG], MK[d] + "-", color=O, ms=3.2, lw=0.8)
    ax[1].plot(XT, [g(n0, f"ratio.R1.d{d}.t{tn[t]}") for t in TG], MK[d] + "--", color=B, mfc="none", ms=3.2, lw=0.8); ax[1].plot(XT, [g(n4, f"bw.ratio.R1.d{d}.t{tn[t]}") for t in TG], MK[d] + "-", color=O, ms=3.2, lw=0.8)
ax[0].axhspan(93, 97, color="0.85", lw=0); ax[0].axhline(90, color="0.4", lw=0.6, ls=":"); ax[0].set_ylim(-3, 103); ax[0].set_ylabel("R1 coverage (%)")
ax[1].axhline(1, color="0.4", lw=0.6, ls=":"); ax[1].set_yscale("log"); ax[1].set_ylabel("R1 prediction/reference")
for a, t in zip(ax, "ab"): a.set_xscale("log"); a.invert_xaxis(); a.text(0.0, 1.03, f"({t})", transform=a.transAxes)
[a.set_xlabel("target $p$") for a in ax]
H = [Line2D([], [], color=B, ls="--", marker="o", mfc="none", ms=3, lw=0.8, label="first window"), Line2D([], [], color=O, ls="-", marker="o", ms=3, lw=0.8, label="second window")] + [Line2D([], [], color="0.2", ls="", marker=MK[d], mfc="none", ms=3.2, label=f"$d={d}$") for d in DS]
fig.legend(handles=H, loc="lower center", ncol=5, frameon=False, fontsize=6, columnspacing=1.0, handletextpad=0.4)
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig("figs/fig4_windows.pdf"); fig.savefig("figs/fig4_windows.png", dpi=200)
