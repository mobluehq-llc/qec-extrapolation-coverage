# CORRECTIVE checks for the Gemini red-team of v5 (PLAN.md "CORRECTIVE AMENDMENT (Gemini red-team v5 checks)").
# Stored data only. Run from this folder (code/):  python3 redteam_checks_gemini.py   (needs numpy, scipy; run on a compute host)
import json, math, csv, sys, glob, os, numpy as np
import coverage as C1
import coverage_bw as C2
from analyze_main import load_truth
import certified
R = '../results'
TG = [1e-3, 1e-4, 1e-5, 1e-6]
out = {}

def eps(P, d): return -0.5 * math.expm1(math.log1p(-2 * P) / d)
def Pinv(e, d): return -0.5 * math.expm1(d * math.log1p(-2 * e))

T = load_truth(); truth = {(d, p): x['truth'] for (d, p), x in T.items()}

# ---------------------------------------------------------------- item 1 (a): reference-value distance factor
a = {}
for tp in TG:
    p3, p5, p7 = truth[(3, tp)], truth[(5, tp)], truth[(7, tp)]
    e3, e5, e7 = eps(p3, 3), eps(p5, 5), eps(p7, 7)
    proj_blk = p5 ** 2 / p3; proj_pr = Pinv(e5 ** 2 / e3, 7)
    a[str(tp)] = dict(p3=p3, p5=p5, p7=p7, eps3=e3, eps5=e5, eps7=e7, lam_blk_35=p3 / p5, lam_blk_57=p5 / p7, lam_pr_35=e3 / e5, lam_pr_57=e5 / e7,
                      factor_block=proj_blk / p7, factor_perround=proj_pr / p7, alg_25_21=25 / 21, factor_block_over_25_21=proj_blk / p7 / (25 / 21))
out['item1a_distance_factor'] = a

# ---------------------------------------------------------------- item 1 (b): R5 per replicate, both windows, shifted interval
def r5_rows(mod, path, label, frozen_csv):
    reps = mod.load_reps(path); res = {}
    hits_blk = {tp: 0 for tp in TG}; hits_pr = {tp: 0 for tp in TG}; rat_blk = {tp: [] for tp in TG}; rat_pr = {tp: [] for tp in TG}; shifts = {tp: [] for tp in TG}
    for rep in reps:
        for tp in TG:
            o = mod.fit_all({d: rep[d] for d in (3, 5, 7)}, tp)
            m3, m5 = o[('R1', 3)][0], o[('R1', 5)][0]; mu, var = o[('R5', 7)]
            P3, P5 = math.exp(m3), math.exp(m5); pt_blk = math.exp(mu)
            assert abs(pt_blk - P5 ** 2 / P3) / pt_blk < 1e-9
            pt_pr = Pinv(eps(P5, 5) ** 2 / eps(P3, 3), 7); f = pt_pr / pt_blk; sd = math.sqrt(var); t = truth[(7, tp)]
            hits_blk[tp] += int(abs(math.log(t) - mu) <= 1.96 * sd)
            hits_pr[tp] += int(abs(math.log(t) - (mu + math.log(f))) <= 1.96 * sd)
            rat_blk[tp].append(pt_blk / t); rat_pr[tp].append(pt_pr / t); shifts[tp].append(f)
    ref = {(r['recipe'], int(r['d']), float(r['target'])): int(r['hits']) for r in csv.DictReader(open(frozen_csv))}
    bad = [(tp, hits_blk[tp], ref[('R5', 7, tp)]) for tp in TG if hits_blk[tp] != ref[('R5', 7, tp)]]
    res['sanity_mismatch_vs_frozen_csv'] = bad
    for tp in TG:
        lo, hi = C1.cp(hits_pr[tp], len(reps)); lo0, hi0 = C1.cp(hits_blk[tp], len(reps))
        res[str(tp)] = dict(reps=len(reps), hits_block=hits_blk[tp], cov_block=hits_blk[tp] / len(reps), hits_perround=hits_pr[tp], cov_perround=hits_pr[tp] / len(reps),
                            cp_perround=[lo, hi], med_ratio_block=float(np.median(rat_blk[tp])), med_ratio_perround=float(np.median(rat_pr[tp])), med_shift=float(np.median(shifts[tp])))
    print(label, 'R5 sanity mismatches:', bad)
    return res
out['item1b_R5_window1'] = r5_rows(C1, f'{R}/fit_d357.json', 'window1', f'{R}/L_coverage.csv')
out['item1b_R5_bw'] = r5_rows(C2, f'{R}/fit_bw.json', 'bw', f'{R}/L_bw_coverage.csv')

# ---------------------------------------------------------------- item 1 (c): d* table, block vs per-round
preds = json.load(open(f'{R}/L_amend_gemini_preds.json')); tp = 1e-3
def med(rc, d): return float(np.median([x[0] for x in preds[f'{rc}|{d}|{tp}|primary']]))
def dstar_blk(p7, lam):
    d = 7
    while p7 * lam ** (-(d - 7) / 2) > 1e-12 and d < 201: d += 2
    return d
def dstar_pr(e7, lam):
    d = 7
    while Pinv(e7 * lam ** (-(d - 7) / 2), d) > 1e-12 and d < 201: d += 2
    return d
frozen = json.load(open(f'{R}/L_amend_gemini.json'))['dstar']
rows = {}
src = {'truth': (truth[(5, tp)], truth[(7, tp)], None)}
for rc in ('R1', 'R2', 'R3', 'R4', 'R4p', 'R2lin'): src[rc] = (med(rc, 5), med(rc, 7), None)
for k, (m5, m7, _) in src.items():
    lb = m5 / m7; db = dstar_blk(m7, lb); lp = eps(m5, 5) / eps(m7, 7); dp = dstar_pr(eps(m7, 7), lp)
    rows[k] = dict(p7=m7, lam_block=lb, dstar_block=db, lam_perround=lp, dstar_perround=dp, frozen_dstar=frozen[k]['dstar'], frozen_lam=frozen[k]['lam'])
# R5: Lambda from R1 d=3,5 medians; anchor = median over replicates of its own d=7 value
m3, m5 = med('R1', 3), med('R1', 5)
p7b = med('R5', 7)
p7p = float(np.median([Pinv(eps(x5[0], 5) ** 2 / eps(x3[0], 3), 7) for x3, x5 in zip(preds[f'R1|3|{tp}|primary'], preds[f'R1|5|{tp}|primary'])]))
lb = m3 / m5; lp = eps(m3, 3) / eps(m5, 5)
rows['R5'] = dict(p7=p7b, lam_block=lb, dstar_block=dstar_blk(p7b, lb), p7_perround=p7p, lam_perround=lp, dstar_perround=dstar_pr(eps(p7p, 7), lp), frozen_dstar=frozen['R5']['dstar'], frozen_lam=frozen['R5']['lam'])
for k, v in rows.items():
    v['dstar_perround_qubits'] = 2 * v['dstar_perround'] ** 2; v['dstar_perround_qubit_rounds'] = 2 * v['dstar_perround'] ** 3
    v['dstar_block_matches_frozen'] = (v['dstar_block'] == v['frozen_dstar'])
out['item1c_dstar'] = rows

# ---------------------------------------------------------------- item 2: d=7 reference combination and shape
names = {1e-3: '0.001', 1e-4: '0.0001', 1e-5: '1e-05', 1e-6: '1e-06'}; tn = {1e-3: '1e-3', 1e-4: '1e-4', 1e-5: '1e-5', 1e-6: '1e-6'}
runs = {}
for p in TG:
    runs[p] = [json.load(open(f'{R}/truth_d7_p{tn[p]}.json')), json.load(open(f'{R}/truthextra_d7_p{names[p]}.json'))]
comb = {}
for p, rs in runs.items():
    t = [x['truth'] for x in rs]; s = [x['se'] for x in rs]; w = [1 / q ** 2 for q in s]
    ivw = sum(a_ * b_ for a_, b_ in zip(w, t)) / sum(w); eq = sum(t) / 2
    assert abs(ivw - truth[(7, p)]) / ivw < 1e-12
    comb[p] = dict(ivw=ivw, ivw_se=sum(w) ** -0.5, eq=eq, eq_se=math.sqrt(sum(q * q for q in s)) / 2, run1=t[0], run2=t[1], se1=s[0], se2=s[1],
                   weight1=w[0] / sum(w), z_run_diff=(t[0] - t[1]) / math.hypot(*s), eq_over_ivw=eq / ivw - 1, secs=[x['secs'] for x in rs])
out['item2a_combination'] = {str(p): v for p, v in comb.items()}
# shape: y_i = p_L/p^4 (units 1e6) for 8 runs ; model A e^-lam (1 + B p)
Y = []; S = []; P = []; L = []
for p, rs in runs.items():
    for x in rs: Y.append(x['truth'] / p ** 4 / 1e6); S.append(x['se'] / p ** 4 / 1e6); P.append(p); L.append(x['lam'])
Y, S, P, L = map(np.array, (Y, S, P, L)); E = np.exp(-L)
def wls(cols, y, s):
    X = np.column_stack(cols); w = 1 / s ** 2; A = X.T @ (X * w[:, None]); Cv = np.linalg.inv(A); b = Cv @ (X.T @ (w * y)); r = (y - X @ b) / s
    return b, Cv, float((r ** 2).sum()), r
from scipy.stats import chi2 as chi2d
b2, C2m, c2, r2 = wls([E, P * E], Y, S)
b1, C1m, c1, r1 = wls([E[2:]], Y[2:], S[2:])        # A only, fitted to p <= 1e-4 (6 values)
pred1_at_1e3 = float(b1[0] * E[:2].mean())          # extrapolated to 1e-3 with B=0
shape = dict(model_A_B=dict(A=float(b2[0]), AB=float(b2[1]), B=float(b2[1] / b2[0]), chi2=c2, dof=6, p_value=float(1 - chi2d.cdf(c2, 6)), z=[float(z) for z in r2]),
             model_A_low=dict(A=float(b1[0]), chi2=c1, dof=5, p_value=float(1 - chi2d.cdf(c1, 5)), z=[float(z) for z in r1]),
             exp_neg_lambda=dict(p1e4_over_p1e6=float(math.exp(-(L[2] - L[6]))), p1e5_over_p1e6=float(math.exp(-(L[4] - L[6])))),
             obs_ratio_ivw=dict(p1e4_over_p1e6=comb[1e-4]['ivw'] / 1e-16 / (comb[1e-6]['ivw'] / 1e-24), p1e5_over_p1e6=comb[1e-5]['ivw'] / 1e-20 / (comb[1e-6]['ivw'] / 1e-24),
                                p1e3_over_p1e6=comb[1e-3]['ivw'] / 1e-12 / (comb[1e-6]['ivw'] / 1e-24)),
             model_A_low_pred_1e3_B0=pred1_at_1e3)
# observed dip significance: ratio 1e-5 (and 1e-4) to the model-expected shape from 1e-6
for q in (1e-4, 1e-5):
    exp = math.exp(-(sum(L[P == q]) / 2 - sum(L[P == 1e-6]) / 2)); obs = (comb[q]['ivw'] / q ** 4) / (comb[1e-6]['ivw'] / 1e-24)
    se = obs * math.hypot(comb[q]['ivw_se'] / comb[q]['ivw'], comb[1e-6]['ivw_se'] / comb[1e-6]['ivw'])
    shape[f'dip_{q:g}_vs_1e-6'] = dict(observed_ratio=obs, structural_expectation_exp_neg_dlam=exp, z=(obs - exp) / se, rel_se_of_ratio=se / obs)
out['item2b_shape'] = shape
smooth = {}
for p in TG:
    lam_p = float(L[P == p][0]); smooth[p] = (float(b1[0] * math.exp(-lam_p)) * 1e6 * p ** 4) if p <= 1e-4 else comb[p]['ivw']   # registered (A,B) model failed (chi2 44/6): post-hoc fallback = A only, fitted to p<=1e-4; p=1e-3 left at its IVW value
out['item2b_smooth_ref'] = {str(p): v for p, v in smooth.items()}
# item 2c: R2 d=7 coverage vs different references
def r2_cov(mod, path, frozen_csv, label):
    reps = mod.load_reps(path); ref = {(r['recipe'], int(r['d']), float(r['target'])): int(r['hits']) for r in csv.DictReader(open(frozen_csv))}
    var = {'ivw': {p: comb[p]['ivw'] for p in TG}, 'equal': {p: comb[p]['eq'] for p in TG}, 'run1': {p: comb[p]['run1'] for p in TG}, 'run2': {p: comb[p]['run2'] for p in TG},
           'smoothed_model': {p: smooth[p] for p in TG}}
    hits = {k: {p: 0 for p in TG} for k in var}; ratio = {k: {p: [] for p in TG} for k in var}
    for rep in reps:
        for tp in TG:
            mu, v = mod.fit_all({d: rep[d] for d in (3, 5, 7)}, tp)[('R2', 7)]
            for k in var:
                hits[k][tp] += int(abs(math.log(var[k][tp]) - mu) <= 1.96 * math.sqrt(v)); ratio[k][tp].append(math.exp(mu) / var[k][tp])
    bad = [(tp, hits['ivw'][tp], ref[('R2', 7, tp)]) for tp in TG if hits['ivw'][tp] != ref[('R2', 7, tp)]]
    print(label, 'R2 d=7 sanity mismatches vs frozen csv:', bad)
    return dict(sanity_mismatch=bad, reps=len(reps), hits={k: {str(p): h for p, h in d_.items()} for k, d_ in hits.items()},
                med_ratio={k: {str(p): float(np.median(x)) for p, x in d_.items()} for k, d_ in ratio.items()})
out['item2c_R2_d7_window1'] = r2_cov(C1, f'{R}/fit_d357.json', f'{R}/L_coverage.csv', 'window1')
out['item2c_R2_d7_bw'] = r2_cov(C2, f'{R}/fit_bw.json', f'{R}/L_bw_coverage.csv', 'bw')

# ---------------------------------------------------------------- item 3: planning numbers (does NOT sample): needs the model's P_w
import truth as TR
plan = {}
for p in TG:
    m = TR.Model(7, p, 1, 'full'); Pw = m.Pw; st = {}
    pooled = {}
    for x in runs[p]:
        for w, (S1, S2, M, nf) in x['strata'].items():
            q = pooled.setdefault(int(w), [0.0, 0.0, 0, 0]); q[0] += S1; q[1] += S2; q[2] += M; q[3] += nf
    a_w = {}; mu = 0.0
    for w, (S1, S2, M, nf) in pooled.items():
        f = S1 / M; s2 = max(S2 / M - f * f, 0.0); mu += Pw[w] * f; a_w[w] = Pw[w] * math.sqrt(s2)
    # w<4 fixed at 0 (as frozen runs; 0 failures): drop. Neyman: M_w ~ a_w; var_total = (sum a_w)^2 / Mtot
    sa = sum(v for w, v in a_w.items() if w >= 4)
    rel = 0.015
    Mtot_plain = (sa / (rel * mu)) ** 2
    # conservative frozen-style se: P_w sqrt((S2+Wmax^2)/M^2) sqrt(max(1,30/ess)); recompute pooled with M-scaling
    Wmax = 1 / (1 - TR.ALPHA); a_c = {}
    for w, (S1, S2, M, nf) in pooled.items():
        if w < 4: continue
        ess = (S1 * S1 / S2) if S2 > 0 else 0.0; a_c[w] = Pw[w] * math.sqrt((S2 / M + Wmax ** 2 / M)) * max(1.0, math.sqrt(30.0 / max(ess, 1.0)))
    sac = sum(a_c.values()); Mtot_cons = (sac / (rel * mu)) ** 2
    thr = 1.25e4  # draws/s/core: totM/secs of the stored runs
    rates = [sum(v[2] for v in x['strata'].values()) / x['secs'] for rs_ in runs.values() for x in rs_]
    plan[str(p)] = dict(mu_pooled=mu, mu_ivw=comb[p]['ivw'], share_by_stratum={str(w): float(Pw[w] * pooled[w][0] / pooled[w][2] / mu) for w in sorted(pooled) if Pw[w] * pooled[w][0] / pooled[w][2] / mu > 0.01},
                        Mtot_plain_1p5pct=Mtot_plain, Mtot_plain_x1p5_safety=1.5 * Mtot_plain, Mtot_conservative_1p5pct=Mtot_cons,
                        top_strata_alloc_fraction={str(w): float(a_w[w] / sa) for w in sorted(a_w, key=lambda w: -a_w[w])[:5] if w >= 4},
                        existing_total_draws=sum(sum(v[2] for v in x['strata'].values()) for x in runs[p]))
rates = [sum(v[2] for v in x['strata'].values()) / x['secs'] for rs_ in runs.values() for x in rs_]
plan['_rate_draws_per_s_per_core'] = dict(min=min(rates), median=float(np.median(rates)), max=max(rates))
out['item3_plan'] = plan

# ---------------------------------------------------------------- item 4: S16
c = json.load(open(f'{R}/cert_d9_result.json')); per = c['A0']['per']; S_ = len(per); dp = 0.05 / (2 * S_)
s16 = {'S': S_, 'dprime': dp, 'sumL': sum(v['P'] * v['lcb'] for v in per.values()), 'sumU_terms': {w: v['P'] * v['ucb'] for w, v in per.items()}, 'tail': c['A0']['tail'], 'L': c['A0']['L'], 'U': c['A0']['U']}
s16['sumU'] = sum(s16['sumU_terms'].values()) + c['A0']['tail']
m9 = certified.truth.Model(9, 1e-4, 1, 'full'); bb = {}
for w in ('11', '12', '10'):
    fs = sorted(glob.glob(f'{R}/cert_d9_chunks/w{w}_c*.npz')); v = np.concatenate([np.load(f)['v'] for f in fs]); n = sum(int(np.load(f)['n']) for f in fs)
    b = certified.bw(m9, int(w)); l, u = certified.bounds(v, n, b, dp); l0, u0 = certified.bounds(np.empty(0), n, b, dp)
    bb[w] = dict(n=n, nfail=len(v), b=b, Ymax=float(v.max()) if len(v) else 0.0, Ymean_of_failing=float(v.mean()) if len(v) else 0.0, f_hat=float(v.sum() / n), ucb_actual=u, ucb_all_zero_data=u0,
                 ucb_stored=per[w]['ucb'], ucb_over_all_zero=u / u0 if u0 else None, P=per[w]['P'], contrib_U=per[w]['P'] * u, lcb=l)
    # Gemini's scale term b ln(2S/delta)/M
    bb[w]['gemini_scale_term'] = b * math.log(1 / dp) / n; bb[w]['grid_lambda_max_term'] = b * (math.log(1 / dp) + math.log(18)) / n
s16['strata'] = bb; out['item4_S16'] = s16

# ---------------------------------------------------------------- item 5: crossings
crs = {}
reps = C1.load_reps(f'{R}/fit_d357.json'); PS = C1.PS
for stat, fn in (('mean', lambda a_: a_.mean(0)), ('median', lambda a_: np.median(a_, 0))):
    Wm = {d: fn(np.array([r[d][0] / r[d][1] for r in reps])) for d in (3, 5, 7)}
    for (x, y) in ((3, 5), (3, 7), (5, 7)):
        ld = np.log(Wm[x]) - np.log(Wm[y]); cr = [math.exp(math.log(PS[i]) + ld[i] / (ld[i] - ld[i + 1]) * (math.log(PS[i + 1]) - math.log(PS[i]))) for i in range(4) if ld[i] * ld[i + 1] < 0]
        crs[f'{stat}_d{x}_d{y}'] = [c_ * 1e3 for c_ in cr]
out['item5_crossings_1e-3'] = crs

# ---------------------------------------------------------------- (b) surrogate: certified bracket around the exact polynomial
sys.path.insert(0, os.path.expanduser('work/L_exact')); sys.path.insert(0, '../exact')
import importlib.util
spec = importlib.util.spec_from_file_location('common', '../exact/common.py'); common = importlib.util.module_from_spec(spec); spec.loader.exec_module(common)
def pw_le4(q):
    c = np.zeros(5); c[0] = 1.0
    for qi in q:
        c[1:] = c[1:] * (1 - qi) + c[:-1] * qi; c[0] *= (1 - qi)
    return c
sb = {}
for d in (3, 5):
    M = common.build(d); m = M['m']; man = json.load(open(f'../exact/d{d}_manifest.json')); mmax = float(m.max())
    for p in (1e-3, 1e-4):
        q = m * p; E = float(np.prod(1 - q)); pw = pw_le4(q); tail5 = float(1 - pw.sum()); lam = float(q.sum())
        Lx = Ux = 0.0; poly4 = sum(v['coef'] * p ** int(k) for k, v in man['coefficients_p_L=sum c_k p^k'].items())
        low = 0.0; extra = 0.0
        for w in range(1, 5):
            Tw = man['shells'][str(w)]['T']; assert man['shells'][str(w)]['complete']
            low += E * p ** w * sum(p ** j * Tw[j] for j in range(4))
            sig = p * w * mmax; assert sig < 1
            extra += E * p ** w * Tw[0] * sig ** 4 / (1 - sig)
        Lx = low; Ux = low + extra + tail5
        # also bracket if the first four shells' own series remainder were the only uncertainty (needs shell>=5 failing sets): reported as P(W>=5)/Lx
        sb[f'd{d}_p{p:g}'] = dict(n_mech=int(len(m)), lam=lam, poly4=poly4, L_exact=Lx, U_exact=Ux, series_extra=extra, P_W_ge5=tail5,
                                   known_higher_order_over_poly4=Lx / poly4 - 1, rel_width=(Ux - Lx) / poly4, rel_width_without_shell5_tail=extra / poly4)
        # re-score the stored certified replicate intervals
        fs = sorted(glob.glob(f'{R}/cert_validation/d{d}_p{p:.0e}'.replace('e-0', 'e-') + f'/cert_exact_d{d}_p{p:g}_r*.json')); rs = [json.load(open(f)) for f in fs]
        sb[f'd{d}_p{p:g}'].update(n_reps=len(rs), inside_poly4=sum(r['L'] <= poly4 <= r['U'] for r in rs),
                                  bracket_inside=sum(r['L'] <= Lx and Ux <= r['U'] for r in rs), bracket_intersects=sum(r['L'] <= Ux and Lx <= r['U'] for r in rs),
                                  lowerbracket_point_inside=sum(r['L'] <= Lx <= r['U'] for r in rs))
        isf = sorted(glob.glob(os.path.expanduser(f'work/L_is/out/d{d}_p{p:g}_s*.json'))); isr = [json.load(open(f)) for f in isf if not f.endswith('.tmp')]
        def dist(lo, hi, a, b): return max(0.0, a - hi, lo - b)
        sb[f'd{d}_p{p:g}'].update(is_runs=len(isr), is_inside_poly4=sum(abs(r['truth'] - poly4) <= 1.96 * r['se'] for r in isr),
                                  is_touch_bracket=sum(dist(Lx, Ux, r['truth'] - 1.96 * r['se'], r['truth'] + 1.96 * r['se']) == 0 for r in isr),
                                  is_touch_lowerpart=sum(dist(Lx, poly4 if poly4 > Lx else Lx, r['truth'] - 1.96 * r['se'], r['truth'] + 1.96 * r['se']) == 0 for r in isr))
out['itemb_surrogate_bracket'] = sb
# ---------------------------------------------------------------- (c) W>120 tail and the betting grids
tl = {}
for d, p in ((9, 1e-4), (7, 1e-3)):
    mm = TR.Model(d, p, 1, 'full') if d == 7 else m9; lam = float(mm.lam); k = 121
    logb = -lam + k * (1 + math.log(lam / k)); tl[f'd{d}_p{p:g}'] = dict(lam=lam, k=k, chernoff_log10_bound_P_W_ge_121=logb / math.log(10),
        note='P(W>=k) <= exp(-lam) (e lam/k)^k, valid for k>lam, Poisson-binomial with mean lam (Chernoff with the Poisson mgf bound)', exact_P_W_gt_12=float(mm.Pw[13:].sum()) if d == 9 else None)
out['itemc_tail120'] = tl
out['itemc_grids'] = dict(GP_multipliers_g=[float(x) for x in certified.GP], GM_multipliers_g=[float(x) for x in certified.GM],
    rule='lambda_j = g_j / b_w ; K+ (lower limit) uses (1+lambda(Y-m)) and skips lambda_j when lambda_j*m>=1 i.e. admissible only for m < b_w/g_j ; K- (upper limit) uses (1+lambda(m-Y)) and skips lambda_j when lambda_j*(b_w-m)>=1; with g<=0.999 every m in [0,b_w) is admissible. Average of the 18 products must reach 1/delta_prime; delta_prime = delta/(2S).',
    admissible_m_over_b_for_Kplus=[float(1 / x) if x > 1 else 1.0 for x in certified.GP])
# ---------------------------------------------------------------- (a) sweep cost
def interp_pl(d):
    pts = {}
    for f_, key in ((f'{R}/fit_d357.json', 'w1'), (f'{R}/fit_bw.json', 'bw')):
        agg = {}
        for x in json.load(open(f_)):
            if x['d'] == d: a_ = agg.setdefault(x['ip'], [0, 0]); a_[0] += x['k']; a_[1] += x['n']
        for ip, (k_, n_) in agg.items(): pts[(C1.PS if key == 'w1' else C2.PS)[ip]] = k_ / n_
    for (dd, p_), x in truth.items():
        if dd == d: pts[p_] = x
    xs = sorted(pts); lx = np.log(xs); ly = np.log([pts[x] for x in xs])
    def f(p_):
        if p_ < xs[0]:
            sl = (ly[1] - ly[0]) / (lx[1] - lx[0]); return math.exp(ly[0] + sl * (math.log(p_) - lx[0]))
        return math.exp(np.interp(math.log(p_), lx, ly))
    return f
rate = {}
n_ = collections = None
import collections as _c
nn = _c.defaultdict(float); ss = _c.defaultdict(float)
for l in open('../logs/L_cov_manifest.jsonl'):
    x = json.loads(l)
    if x['stage'] == 'S1': nn[x['res']['d']] += x['res']['n']; ss[x['res']['d']] += x['secs']
rate = {d: nn[d] / ss[d] for d in nn}
sw = {}; RATIOS = [0.2, 0.3, 0.45, 0.7, 1.0, 1.5, 2.0]; SPAN = 2.67; NREP = 100
for d, px in ((3, 7.74e-3), (5, 7.74e-3)):
    f = interp_pl(d); rows_ = {}
    for r in RATIOS:
        pmax = r * px; pts = [pmax / SPAN ** (i / 4) for i in range(5)]
        shots = sum(300 / f(p_) for p_ in pts); rows_[str(r)] = dict(pmax=pmax, shots_per_rep=shots, cpu_h_100reps_at_measured_rate=shots * NREP / rate[d] / 3600,
                                                                      cpu_h_100reps_at_third_of_rate=shots * NREP / (rate[d] / 3) / 3600)
    sw[f'd{d}'] = dict(rate_shots_per_s=rate[d], windows=rows_, total_cpu_h_measured=sum(v['cpu_h_100reps_at_measured_rate'] for v in rows_.values()),
                       total_cpu_h_third=sum(v['cpu_h_100reps_at_third_of_rate'] for v in rows_.values()))
out['itema_sweep_cost'] = sw
json.dump(out, open(f'{R}/L_redteam_gemini.json', 'w'), indent=1, default=float)
print('wrote', f'{R}/L_redteam_gemini.json')
