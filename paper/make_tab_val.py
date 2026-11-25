"""tab_v9_val.tex: the 16 direct-Monte-Carlo validation cells of the reference estimator with their z-scores (supplement, arXiv v2).
Reads results/L_truth_validation.csv only. Run from this folder: python3 make_tab_val.py"""
import csv, math
R = list(csv.DictReader(open("../results/L_truth_validation.csv")))
def sci(x, nd=3):
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    if round(m, nd - 1) >= 10: m /= 10; e += 1
    return f"{m:.{nd-1}f}\\times10^{{{e}}}"
assert len(R) == 16 and all(r["eligible"] == "1" and r["ok"] == "1" for r in R) and all(abs(float(r["z"])) <= 2 for r in R)
with open("tab_v9_val.tex", "w") as f:
    f.write("\\begin{tabular}{rrrrrrr}\\hline\\hline\n$d$ & $p$ & estimator & direct Monte Carlo & logical errors & shots & $z$\\\\\\hline\n")
    for r in sorted(R, key=lambda r: (int(r["d"]), float(r["p"]))):
        f.write(f"{r['d']} & ${sci(float(r['p']), 2)}$ & ${sci(float(r['truth']))}$ & ${sci(float(r['direct']))}$ & {r['direct_errors']} & ${sci(float(r['shots']), 2)}$ & ${float(r['z']):+.2f}$\\\\\n")
    f.write("\\hline\\hline\n\\end{tabular}\n")
print("max |z|", max(abs(float(r["z"])) for r in R))
