# Exact reference-range extrema (ChatGPT final round, 2026-11-26). Replaces the 81-point grid of amend_gemini.py / r4p_bw.py.
# The hit count k(t) = #{i: lo_i <= t <= hi_i} is piecewise constant in t, changing only at interval endpoints. Its max over
# [t0-1.96se, t0+1.96se] is attained at a range end or an endpoint inside the range; its min at a range end or at a midpoint between
# consecutive breakpoints. Evaluating all of these gives the exact extrema. First window: stored preds (L_amend_gemini_preds.json).
# Second window: preds rebuilt by r4p_bw.py's own fits (no simulation). Writes results/L_exact_ref_range.json. Run from code/.
import json, numpy as np
import r4p_bw as B
from coverage_bw import load_reps, fit_all, TARGETS
from analyze_bw import load_truth as lt_bw
from analyze_main import load_truth as lt_w1
R = '../results'
def ext(A, t0, se):
    a, b = t0 - 1.96 * se, t0 + 1.96 * se
    br = np.unique(np.r_[a, b, A[:, 1][(A[:, 1] >= a) & (A[:, 1] <= b)], A[:, 2][(A[:, 2] >= a) & (A[:, 2] <= b)]])
    cand = np.r_[br, (br[:-1] + br[1:]) / 2]
    k = np.array([int(((A[:, 1] <= t) & (t <= A[:, 2])).sum()) for t in cand])
    g = np.linspace(a, b, 81); kg = np.array([int(((A[:, 1] <= t) & (t <= A[:, 2])).sum()) for t in g])
    return int(k.min()), int(k.max()), int(kg.min()), int(kg.max())
out = {}
T1 = lt_w1(); P1 = json.load(open(f'{R}/L_amend_gemini_preds.json'))
for key, L in P1.items():
    rc, d, tp, sc = key.split('|'); x = T1[(int(d), float(tp))]
    out['w1|' + key] = dict(zip(('kmin', 'kmax', 'kmin_grid', 'kmax_grid'), ext(np.array(L), x['truth'], x['se'])))
T2 = lt_bw(); reps = load_reps(f'{R}/fit_bw.json'); pred = {}
for rep in reps:
    for tp in TARGETS:
        for scaled in (False, True):
            dd = {d: rep[d] for d in B.DS}
            for (rc, d), (mu, var) in fit_all(dd, tp, scaled=scaled).items():
                sd = np.sqrt(var); pred.setdefault((rc, d, tp, scaled), []).append((np.exp(mu), np.exp(mu - 1.96 * sd), np.exp(mu + 1.96 * sd)))
            for (rc, d), v in B.fit_amend(dd, tp, scaled=scaled).items(): pred.setdefault((rc, d, tp, scaled), []).append(v)
for (rc, d, tp, sc), L in pred.items():
    x = T2[(d, tp)]
    out[f'w2|{rc}|{d}|{tp}|{"scaled" if sc else "primary"}'] = dict(zip(('kmin', 'kmax', 'kmin_grid', 'kmax_grid'), ext(np.array(L), x['truth'], x['se'])))
json.dump(out, open(f'{R}/L_exact_ref_range.json', 'w'), indent=1)
ch = {k: v for k, v in out.items() if (v['kmin'], v['kmax']) != (v['kmin_grid'], v['kmax_grid'])}
print('cells', len(out), 'changed', len(ch)); [print(k, v) for k, v in sorted(ch.items())]
