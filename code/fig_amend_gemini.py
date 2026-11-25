# Plotting only (no inference). Figures for the Gemini corrective amendment, from stored files. Run from code/: <venv python> fig_amend_gemini.py
import json, glob, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from coverage import fit_all, load_reps, PS
from analyze_main import load_truth
R = '../results'; O = '../figs_amend'; os.makedirs(O, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.spines.top': False, 'axes.spines.right': False})
T = load_truth(); reps = load_reps(f'{R}/fit_d357.json'); DS = (3, 5, 7)
W = {d: np.array([r[d][0] / r[d][1] for r in reps]) for d in DS}      # reps x 5
Wm = {d: np.median(W[d], 0) for d in DS}
ex = {3: {2: 474.63116381773466, 3: -8719.172035221301, 4: -35821.32144124294}, 5: {3: 62996.22845669943, 4: -684898.4807474464}}
col = {3: '#0072B2', 5: '#D55E00', 7: '#009E73'}
# R2 fitted prefactor (median over replicates of exp(b0) at the registered fit)
pref = {}
for d in DS:
    v = [np.exp(mu - ((d + 1) / 2) * np.log(1e-3)) for rep in reps for (rc, dd), (mu, var) in fit_all({x: rep[x] for x in DS}, 1e-3).items() if rc == 'R2' and dd == d]
    pref[d] = float(np.median(v))
# finite-size crossing: d=5 vs d=7 window medians (log-linear in ln p)
ld = np.log(Wm[5]) - np.log(Wm[7]); cross = None
for i in range(4):
    if ld[i] * ld[i + 1] < 0:
        f = ld[i] / (ld[i] - ld[i + 1]); cross = float(np.exp(np.log(PS[i]) + f * (np.log(PS[i + 1]) - np.log(PS[i]))))
print('prefactors (R2, median over reps)', pref, ' d5/d7 window-median crossing', cross, ' ln ratio per window point', ld.tolist())
# ---- Fig A
fig, ax = plt.subplots(figsize=(6.2, 3.6)); pg = np.logspace(-6, np.log10(1.5e-2), 300)
for d in DS:
    h = (d + 1) / 2; tp = sorted((p, x['truth'], x['se']) for (dd, p), x in T.items() if dd == d)
    ax.errorbar([a[0] for a in tp], [a[1] / a[0] ** h for a in tp], yerr=[1.96 * a[2] / a[0] ** h for a in tp], fmt='o', ms=3.5, color=col[d], label=f'd={d} reference estimates (±1.96 se)')
    ax.plot(PS, Wm[d] / PS ** h, 's', mfc='none', ms=6, color=col[d], label='window medians (300 replicates)' if d == 3 else None)
    ax.axhline(pref[d], color=col[d], ls=':', lw=1)
    if d in ex:
        y = sum(c * pg ** (k - h) for k, c in ex[d].items()); ok = len(ex[d]) * 0 + 1
        ax.plot(pg, y, '-', color=col[d], lw=1, alpha=.7, label='exact polynomial of the surrogate ($d=3,5$)' if d == 3 else None)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(8e-7, 1.8e-2); ax.set_xlabel('physical error rate p'); ax.set_ylabel(r'$p_L / p^{(d+1)/2}$')
ax.set_ylim(1e2, 2e9)  # caption carries the notes (dotted = R2 fitted prefactor; surrogate polynomials truncated; d=7 not drawn)
ax.legend(fontsize=6, loc='upper left', ncol=2, frameon=False); fig.savefig(f'{O}/figA_prefactor.png', dpi=220, bbox_inches='tight'); fig.savefig(f'{O}/figA_prefactor.pdf', bbox_inches='tight'); plt.close(fig)
# ---- Fig B money plot
fig, axs = plt.subplots(1, 3, figsize=(9, 3.4), sharey=True)
lo = {(rc, d): [] for rc in ('R1', 'R2') for d in DS}; hi = {k: [] for k in lo}; mid = {k: [] for k in lo}
for tp in pg:
    acc = {k: [] for k in lo}
    for rep in reps:
        for (rc, d), (mu, var) in fit_all({x: rep[x] for x in DS}, tp).items():
            if rc in ('R1', 'R2'): acc[(rc, d)].append((mu - 1.96 * np.sqrt(var), mu, mu + 1.96 * np.sqrt(var)))
    for k, v in acc.items(): a = np.array(v); lo[k].append(np.median(a[:, 0])); mid[k].append(np.median(a[:, 1])); hi[k].append(np.median(a[:, 2]))
for ax, d in zip(axs, DS):
    for rc, c in (('R1', '#CC79A7'), ('R2', '#0072B2')):
        ax.fill_between(pg, np.exp(lo[(rc, d)]), np.exp(hi[(rc, d)]), color=c, alpha=.35, lw=0, label=f'{rc} median 95% interval')
        ax.plot(pg, np.exp(mid[(rc, d)]), color=c, lw=.8)
    tp = sorted((p, x['truth'], x['se']) for (dd, p), x in T.items() if dd == d)
    ax.errorbar([a[0] for a in tp], [a[1] for a in tp], yerr=[1.96 * a[2] for a in tp], fmt='ko', ms=3, label='reference estimates')
    ax.plot(PS, Wm[d], 's', mfc='none', color='r', ms=6, label='window medians')
    if cross: ax.axvline(cross, color='gray', ls='--', lw=.8, label=f'finite-size crossing d5/d7 ≈ {cross:.1e}')
    ax.axvspan(4e-3, 1.5e-2, color='gray', alpha=.08); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-19, 1); ax.set_xlim(8e-7, 1.8e-2)
    ax.set_title(f'd={d}'); ax.set_xlabel('p')
axs[0].set_ylabel(r'$p_L$'); axs[0].legend(fontsize=5.5, loc='lower right'); fig.savefig(f'{O}/figB_money.png', dpi=220, bbox_inches='tight'); fig.savefig(f'{O}/figB_money.pdf', bbox_inches='tight')
