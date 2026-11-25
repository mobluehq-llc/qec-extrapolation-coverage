# 2026-11-26: write the exact reference-range extrema (results/L_exact_ref_range.json, code/exact_ref_range.py) into the stored
# L_amend_gemini.json (first window) and L_r4p_bw.json (second window); the 81-point-grid values are kept as hit_min_grid81 / hit_max_grid81.
import json
R = '../results'; E = json.load(open(f'{R}/L_exact_ref_range.json'))
for w, f, sub in (('w1', 'L_amend_gemini.json', 'reference_uncertainty'), ('w2', 'L_r4p_bw.json', None)):
    J = json.load(open(f'{R}/{f}')); D = J[sub] if sub else J; n = 0
    for k, v in D.items():
        if k.startswith('_') or 'hit_min_grid81' in v: continue
        e = E[f'{w}|{k}']; assert (round(v['hit_min'] * 300), round(v['hit_max'] * 300)) == (e['kmin_grid'], e['kmax_grid']), k
        v['hit_min_grid81'], v['hit_max_grid81'] = v['hit_min'], v['hit_max']; v['hit_min'], v['hit_max'] = e['kmin'] / 300, e['kmax'] / 300; n += (e['kmin'], e['kmax']) != (e['kmin_grid'], e['kmax_grid'])
    json.dump(J, open(f'{R}/{f}', 'w'), indent=1); print(f, 'changed cells', n)
