# CORRECTIVE amendment (Gemini review, 2026-11-25): R4' (exponent (d+1)/2), R2-lin, single-point Lambda projection,
# d* table. Uses ONLY stored window data (results/fit_d357.json) and stored truth cells. No simulation.
# Run from this folder:  python3 amend_gemini.py   (needs numpy, scipy)
import json, sys, numpy as np
from coverage import wls, load_reps, cp, fit_all, PS, X0, TARGETS
from analyze_main import load_truth
R = '../results'
DS = (3, 5, 7)
def fit_amend(data, tp, scaled=False):
    """returns {(recipe,d):(center, halfwidth_kind, ...)} -> dict recipe,d -> (pred, lo, hi) with lo/hi the 95% interval on the p_L scale"""
    out = {}
    xs = np.log(tp)
    def sc(C, chi2, dof): return C * max(1.0, chi2 / dof) if scaled and dof > 0 else C
    for d in data:
        k, n = data[d]; P = k / n; y = np.log(P); w = k / (1 - P); h = (d + 1) / 2
        # R4' : p^h * exp(c0 + c1 p + c2 p^2)   (same WLS-on-ln form as registered R4, exponent corrected)
        X = np.c_[np.ones(5), PS / 1e-2, (PS / 1e-2) ** 2]
        b, C, c2, dof = wls(X, y - h * X0, w); C = sc(C, c2, dof)
        g = np.array([1, tp / 1e-2, (tp / 1e-2) ** 2]); mu = g @ b + h * xs; sd = np.sqrt(g @ C @ g)
        out[('R4p', d)] = (np.exp(mu), np.exp(mu - 1.96 * sd), np.exp(mu + 1.96 * sd))
        # R2-lin : p^h (c0 + c1 p) fitted on the LINEAR scale to z = P/p^h, weights 1/var(z) = k/(z^2 (1-P)) from the observed z
        z = P / PS ** h; s = 1e-2 ** h; zz = z / s; wz = k / (zz ** 2 * (1 - P))
        X = np.c_[np.ones(5), PS / 1e-2]
        b, C, c2, dof = wls(X, zz, wz); C = sc(C, c2, dof)
        g = np.array([1, tp / 1e-2]); pred = (g @ b) * s * tp ** h; sd = np.sqrt(g @ C @ g) * s * tp ** h
        out[('R2lin', d)] = (pred, pred - 1.96 * sd, pred + 1.96 * sd)
    return out
def main():
    T = load_truth(); truth = {(d, p): x['truth'] for (d, p), x in T.items()}
    reps = load_reps(f'{R}/fit_d357.json'); print('replicates', len(reps))
    res = {'coverage': {}, 'ratio': {}, 'lambda_proj': {}, 'dstar': {}}
    # per-rep predictions
    pred = {}   # (recipe,d,tp,scaled)-> list of (pred,lo,hi)
    for rep in reps:
        for tp in TARGETS:
            for scaled in (False, True):
                o = fit_all({d: rep[d] for d in DS}, tp, scaled=scaled)
                for (rc, d), (mu, var) in o.items():
                    sd = np.sqrt(var)
                    pred.setdefault((rc, d, tp, scaled), []).append((np.exp(mu), np.exp(mu - 1.96 * sd), np.exp(mu + 1.96 * sd)))
                a = fit_amend({d: rep[d] for d in DS}, tp, scaled=scaled)
                for (rc, d), v in a.items(): pred.setdefault((rc, d, tp, scaled), []).append(v)
    for (rc, d, tp, scaled), L in sorted(pred.items()):
        t = truth[(d, tp)]; A = np.array(L); hit = ((A[:, 1] <= t) & (t <= A[:, 2])).mean() * 1.0
        k = int(round(hit * len(L))); lo, hi = cp(k, len(L)); r = A[:, 0] / t
        res['coverage'][f'{rc}|{d}|{tp}|{"scaled" if scaled else "primary"}'] = dict(hits=k, reps=len(L), cov=k / len(L), cp_lo=lo, cp_hi=hi,
            med_ratio=float(np.median(r)), q25=float(np.quantile(r, .25)), q75=float(np.quantile(r, .75)), med_width_ratio=float(np.median(A[:, 2] / A[:, 1])) if (A[:, 1] > 0).all() else None)
    # sanity: the re-derived registered R4/R2 coverage must equal the frozen CSV (R4: coverage.py; here via fit_all intervals on p_L scale == ln-scale test)
    import csv
    ref = {(r['recipe'], int(r['d']), float(r['target'])): int(r['hits']) for r in csv.DictReader(open(f'{R}/L_coverage.csv'))}
    bad = [(rc, d, tp, res['coverage'][f'{rc}|{d}|{tp}|primary']['hits'], ref[(rc, d, tp)]) for (rc, d, tp) in ref
           if rc in ('R1', 'R2', 'R3', 'R4', 'R5') and res['coverage'][f'{rc}|{d}|{tp}|primary']['hits'] != ref[(rc, d, tp)]]
    print('SANITY registered-recipe hit mismatches vs frozen L_coverage.csv:', len(bad), bad[:5]); res['sanity_mismatch'] = len(bad)
    # exponent test: R4/R4' log-share
    sh = {}
    for d in DS:
        for tp in TARGETS:
            r4 = np.array([x[0] for x in pred[('R4', d, tp, False)]]); r4p = np.array([x[0] for x in pred[('R4p', d, tp, False)]]); t = truth[(d, tp)]
            m4, m4p = np.median(r4 / t), np.median(r4p / t)
            sh[f'{d}|{tp}'] = dict(R4_over_truth=float(m4), R4p_over_truth=float(m4p), R4_over_R4p=float(np.median(r4 / r4p)),
                                   log_share=float(np.log(np.median(r4 / r4p)) / np.log(m4)) if abs(np.log(m4)) > 1e-9 else None)
    res['ratio'] = sh
    # single operating point Lambda projection at p=4e-3 (ip=0): P7hat = P3 * Lambda35^{-2} = P5^2/P3 ; ln-interval, var=(1-P5)/k5*4+(1-P3)/k3
    t7 = truth[(7, 4e-3)]; hits = 0; rat = []; wid = []
    for rep in reps:
        k3, n3 = rep[3][0][0], rep[3][1][0]; k5, n5 = rep[5][0][0], rep[5][1][0]; P3, P5 = k3 / n3, k5 / n5
        mu = 2 * np.log(P5) - np.log(P3); sd = np.sqrt(4 * (1 - P5) / k5 + (1 - P3) / k3)
        hits += int(abs(np.log(t7) - mu) <= 1.96 * sd); rat.append(np.exp(mu) / t7); wid.append(np.exp(2 * 1.96 * sd))
    lo, hi = cp(hits, len(reps))
    res['lambda_proj'] = dict(p=4e-3, truth_d7=t7, hits=hits, reps=len(reps), cov=hits / len(reps), cp_lo=lo, cp_hi=hi, med_ratio=float(np.median(rat)),
                              q25=float(np.quantile(rat, .25)), q75=float(np.quantile(rat, .75)), med_interval_ratio=float(np.median(wid)))
    # (a) R5 decomposition: reference-rate distance projection p_L(5)^2/p_L(3) vs d=7 reference; then R5 from fitted rates vs that projection
    dec = {}
    for tp in TARGETS:
        proj = truth[(5, tp)] ** 2 / truth[(3, tp)]; r5 = float(np.median([x[0] for x in pred[('R5', 7, tp, False)]]))
        dec[str(tp)] = dict(truth3=truth[(3, tp)], truth5=truth[(5, tp)], truth7=truth[(7, tp)], ref_projection=proj, distance_factor=proj / truth[(7, tp)],
                            R5_median=r5, R5_over_truth=r5 / truth[(7, tp)], fit_input_factor=r5 / proj)
    res['R5_decomposition'] = dec
    # (b) R4 mechanism: successive-decade ratios of median R4/truth and R4/R4' per d, vs sqrt(10)
    mech = {}
    for d in DS:
        a = [sh[f'{d}|{tp}']['R4_over_truth'] for tp in TARGETS]; b = [sh[f'{d}|{tp}']['R4_over_R4p'] for tp in TARGETS]
        mech[str(d)] = dict(R4_over_truth_steps=[a[i + 1] / a[i] for i in range(3)], R4_over_R4p_steps=[b[i + 1] / b[i] for i in range(3)], sqrt10=float(np.sqrt(10)))
    res['R4_mechanism'] = mech
    # (c) reference-uncertainty bounds: ONE shared reference value t varies over truth +/- 1.96 se (linear); min/max hit fraction over an 81-point grid
    def cls(f): return 'miscal' if f < 0.90 else ('calibrated' if 0.93 <= f <= 0.97 else 'other')
    ru = {}
    for (rc, d, tp, scaled), L in sorted(pred.items()):
        A = np.array(L); t0 = truth[(d, tp)]; se = T[(d, tp)]['se']; lo_, hi_ = t0 - 1.96 * se, t0 + 1.96 * se; br = np.unique(np.r_[lo_, hi_, A[:, 1:3].ravel()[(A[:, 1:3].ravel() >= lo_) & (A[:, 1:3].ravel() <= hi_)]]); grid = np.r_[br, (br[:-1] + br[1:]) / 2]  # exact extrema (2026-11-26): every breakpoint and every gap midpoint; was an 81-point grid
        fr = np.array([((A[:, 1] <= t) & (t <= A[:, 2])).mean() for t in grid]); c = {cls(f) for f in fr}
        ru[f'{rc}|{d}|{tp}|{"scaled" if scaled else "primary"}'] = dict(hit_at_point=float(((A[:, 1] <= t0) & (t0 <= A[:, 2])).mean()), hit_min=float(fr.min()), hit_max=float(fr.max()),
                                                                          verdict_classes=sorted(c), can_flip=len(c) > 1, truth_rel_se=se / t0)
    res['reference_uncertainty'] = ru
    # d* table at p=1e-3, target 1e-12, qubits 2 d^2. Extrapolation rule: Lambda_{5->7} = pred5/pred7 (medians over reps) from the recipe's own d=5,7 predictions;
    # R5 (d=7 only): Lambda from its own d=3,5 basis = median R1 d=3 / median R1 d=5. Truth: real-circuit truth at 1e-3 for d=5,7.
    def dstar(p7, lam):
        d = 7
        while p7 * lam ** (-(d - 7) / 2) > 1e-12 and d < 200: d += 2
        return d
    tbl = {}
    tp = 1e-3
    lamT = truth[(5, tp)] / truth[(7, tp)]; tbl['truth'] = dict(p7=truth[(7, tp)], lam=lamT)
    for rc in ('R1', 'R2', 'R3', 'R4', 'R4p', 'R2lin'):
        m = {d: float(np.median([x[0] for x in pred[(rc, d, tp, False)]])) for d in (5, 7)}
        tbl[rc] = dict(p7=m[7], lam=m[5] / m[7])
    m3 = np.median([x[0] for x in pred[('R1', 3, tp, False)]]); m5 = np.median([x[0] for x in pred[('R1', 5, tp, False)]])
    tbl['R5'] = dict(p7=float(np.median([x[0] for x in pred[('R5', 7, tp, False)]])), lam=float(m3 / m5))
    for rc, v in tbl.items():
        v['dstar'] = dstar(v['p7'], v['lam']) if v['lam'] > 1 else None; v['qubits'] = 2 * v['dstar'] ** 2 if v['dstar'] else None; v['rounds'] = v['dstar']; v['spacetime_qubit_rounds'] = 2 * v['dstar'] ** 3 if v['dstar'] else None
    res['dstar'] = tbl
    json.dump(res, open(f'{R}/L_amend_gemini.json', 'w'), indent=1, default=float)
    json.dump({f'{k[0]}|{k[1]}|{k[2]}|{"scaled" if k[3] else "primary"}': [list(map(float, x)) for x in v] for k, v in pred.items()}, open(f'{R}/L_amend_gemini_preds.json', 'w'))
if __name__ == '__main__': main()
