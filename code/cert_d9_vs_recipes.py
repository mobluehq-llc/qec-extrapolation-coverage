"""Pre-registered comparison: certified d=9 interval vs each recipe's 95% interval over the 300 holdout-window replicates."""
import sys, json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from coverage import *
R = os.path.join(HERE, "..", "results"); res = json.load(open(f"{R}/cert_d9_result.json"))
reps3 = load_reps(f"{R}/fit_d357.json"); reps9 = load_reps(f"{R}/fit_d9.json"); reps = [{**a, **b} for a, b in zip(reps3, reps9)]
out = {}
for tag in ("A0", "assumption_free"):
    L, U = res[tag]["L"], res[tag]["U"]; row = {}
    for rc in ("R1", "R2", "R3", "R4", "R5"):
        above = below = 0; lows = []; highs = []
        for data in reps:
            mu, var = fit_all({d: data[d] for d in (3, 5, 9)}, 1e-4)[(rc, 9)]; lo, hi = np.exp(mu - 1.96 * var ** .5), np.exp(mu + 1.96 * var ** .5)
            lows.append(lo); highs.append(hi); above += (U < lo); below += (L > hi)
        row[rc] = dict(certified_missed_high=int(above), certified_missed_low=int(below), n=len(reps), median_lo=float(np.median(lows)), median_hi=float(np.median(highs)))
    out[tag] = dict(L=L, U=U, recipes=row)
    print(tag, "L=%.3g U=%.3g" % (L, U)); [print(" ", rc, v) for rc, v in row.items()]
d = out["A0"]["recipes"]; dec = all(d[r]["certified_missed_high"] + d[r]["certified_missed_low"] >= 285 for r in ("R1", "R5"))
out["decisive_bar_A0_R1_R5_ge_285of300"] = bool(dec); print("DECISIVE (A-0, R1 and R5 >=285/300):", dec)
json.dump(out, open(f"{R}/cert_d9_vs_recipes.json", "w"), indent=1)
