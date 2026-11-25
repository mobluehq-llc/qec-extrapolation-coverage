# L item 8(b): zero-fault baseline E(p) divided out before the R2/R3 fit (RESULT, 2026-11-24)
Pre-registered in PLAN.md item (b); clarification amendment (commit dedb4bb) written before any R2E/R3E number existed. Stored data only, no new simulation. Code `code/zerofault_r2e.py`, log `logs/zerofault_r2e_20261124.txt`, outputs `results/L_zerofault.json`, `results/L_zerofault_coverage.csv`. [internal operations note redacted] Sanity gate passed: plain R2/R3 hit counts re-derived equal the frozen `L_coverage.csv` / `L_bw_coverage.csv` in every cell. E(p) = exp(sum ln(1-q_i)) over undecomposed DEM mechanisms (logE d=3, p=1e-3: -0.1473, equals -lambda to 3 digits; d=5: -0.744).

## Verdicts by the registered bar (all three NOT MET; the prediction was wrong)
- **P1 NOT MET.** First window, d=3,5, R2E+R3E, 16 cells: median prediction/reference in [0.90,1.10] in 3 of 16 (all R3E d=3: 0.96, 0.95, 0.96; its 1e-3 cell is 0.83). Range 0.83 to 5011. R2E: d=3 2.6-3.0, d=5 257-528 (R2 alone: 0.77-0.85).
- **P2 NOT MET.** Point-reference coverage below 90% in 16 of 16 cells (limit 8). Best cell 172/300 (R3E d=3, 1e-4), the rest 0-171.
- **P3 NOT MET.** Second window, R2E d=3,5: robust-below in 8 of 8 cells (limit 2; R2: 3 of 8). R2E hits 0/300 in every cell; median ratio d=3 1.22-1.41, d=5 3.3-6.7.
- **P4 (d=7, reported, not scored):** R2E over-predicts by 28-275x (second window) and 4e6-1e8x (first window); coverage 0/300 everywhere.

## Second window, R3E (reported, not scored)
d=3 ratios 0.90-1.05, hits 2 (1e-3), 160, 235, 161 of 300 (labels: below, below, undecided, below); d=5 ratios 5.9-12.2, 0/300. Chi2-scaled sensitivity (in the CSV, `scaled=1`): R3E d=3 covers 300/300 in all four cells (intervals widen enough), d=5 stays 0/300.

## Reading
Dividing by the exact E(p) does not fix the known part of the drift; it overcorrects. The window data divided by E(p) are P/E = sum_w lambda^w f_w / w! (roughly), which still carries all the higher-weight terms, and those are large in the windows (E at d=5 is 0.05 to 0.0002 in the straddling window), so the fitted prefactor is inflated and the extrapolation, where E is near 1, lands far above the reference. The interval half-width stays about 5% (ln scale), so the miss is certified. Only R3E at d=3 is near unity, and not reliably calibrated (hits 113-235 of 300 at 1e-4 to 1e-6). As pre-registered, this is not a recommendation: the recipe is dropped as a candidate fix. No frozen verdict or number in the paper changes; the .tex was not edited.
