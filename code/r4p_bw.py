# Amendment R4P-BW: R4' and R2-lin on the stored below-threshold window (results/fit_bw.json). No simulation.
import json, csv, numpy as np
from coverage_bw import wls, load_reps, cp, fit_all, PS, X0, TARGETS
from analyze_bw import load_truth
R = '../results'; DS = (3, 5, 7)
def fit_amend(data, tp, scaled=False):
    out = {}; xs = np.log(tp)
    def sc(C, chi2, dof): return C * max(1.0, chi2 / dof) if scaled and dof > 0 else C
    for d in data:
        k, n = data[d]; P = k / n; y = np.log(P); w = k / (1 - P); h = (d + 1) / 2
        X = np.c_[np.ones(5), PS / 1e-2, (PS / 1e-2) ** 2]
        b, C, c2, dof = wls(X, y - h * X0, w); C = sc(C, c2, dof)
        g = np.array([1, tp / 1e-2, (tp / 1e-2) ** 2]); mu = g @ b + h * xs; sd = np.sqrt(g @ C @ g)
        out[('R4p', d)] = (np.exp(mu), np.exp(mu - 1.96 * sd), np.exp(mu + 1.96 * sd))
        z = P / PS ** h; s = 1e-2 ** h; zz = z / s; wz = k / (zz ** 2 * (1 - P))
        X = np.c_[np.ones(5), PS / 1e-2]
        b, C, c2, dof = wls(X, zz, wz); C = sc(C, c2, dof)
        g = np.array([1, tp / 1e-2]); pred = (g @ b) * s * tp ** h; sd = np.sqrt(g @ C @ g) * s * tp ** h
        out[('R2lin', d)] = (pred, pred - 1.96 * sd, pred + 1.96 * sd)
    return out
def band(f): return 'below90' if f < 0.90 else ('inband' if 0.93 <= f <= 0.97 else ('above97' if f > 0.97 else 'gap'))
def main():
    T = load_truth(); truth = {(d, p): x['truth'] for (d, p), x in T.items()}
    reps = load_reps(f'{R}/fit_bw.json'); print('replicates', len(reps))
    pred = {}
    for rep in reps:
        for tp in TARGETS:
            for scaled in (False, True):
                dd = {d: rep[d] for d in DS}
                for (rc, d), (mu, var) in fit_all(dd, tp, scaled=scaled).items():
                    sd = np.sqrt(var); pred.setdefault((rc, d, tp, scaled), []).append((np.exp(mu), np.exp(mu - 1.96 * sd), np.exp(mu + 1.96 * sd)))
                for (rc, d), v in fit_amend(dd, tp, scaled=scaled).items(): pred.setdefault((rc, d, tp, scaled), []).append(v)
    ref = {(r['recipe'], int(r['d']), float(r['target'])): int(r['hits']) for r in csv.DictReader(open(f'{R}/L_bw_coverage.csv'))}
    bad = []
    res = {}
    for (rc, d, tp, scaled), L in sorted(pred.items()):
        A = np.array(L); t = truth[(d, tp)]; hit = ((A[:, 1] <= t) & (t <= A[:, 2])); k = int(hit.sum()); lo, hi = cp(k, len(L))
        if not scaled and rc in ('R1', 'R2', 'R3', 'R4', 'R5') and ref[(rc, d, tp)] != k: bad.append((rc, d, tp, k, ref[(rc, d, tp)]))
        se = T[(d, tp)]['se']; grid = np.linspace(t - 1.96 * se, t + 1.96 * se, 81)
        fr = np.array([((A[:, 1] <= g) & (g <= A[:, 2])).mean() for g in grid]); cl = sorted({band(f) for f in fr})
        res[f'{rc}|{d}|{tp}|{"scaled" if scaled else "primary"}'] = dict(hits=k, reps=len(L), cov=k / len(L), cp_lo=lo, cp_hi=hi, band=band(k / len(L)),
            med_ratio=float(np.median(A[:, 0] / t)), hit_min=float(fr.min()), hit_max=float(fr.max()), classes=cl, truth_rel_se=se / t)
    print('SANITY mismatches vs frozen L_bw_coverage.csv:', len(bad), bad[:5]); res['_sanity_mismatch'] = len(bad)
    json.dump(res, open(f'{R}/L_r4p_bw.json', 'w'), indent=1)
    for rc in ('R4p', 'R2lin', 'R4'):
        for sc in ('primary', 'scaled'):
            for grp, ds in (('d3,5', (3, 5)), ('d7', (7,))):
                c = [res[f'{rc}|{d}|{tp}|{sc}'] for d in ds for tp in TARGETS]
                print(rc, sc, grp, 'inband', sum(x['band'] == 'inband' for x in c), '/', len(c), 'below90', sum(x['band'] == 'below90' for x in c), 'above97', sum(x['band'] == 'above97' for x in c))
if __name__ == '__main__': main()
