# RESULTS: CORRECTIVE amendment (Gemini exponent finding + ChatGPT items), 2026-11-25

Pre-registered in `PLAN.md` ("CORRECTIVE AMENDMENT" commit 36679e0 and "Amendment 1b" commit 70a18c5), both committed BEFORE any number below was computed. Code `code/amend_gemini.py` (sha256 prefix 9e07b882659eba4d, as registered), plotting `code/fig_amend_gemini.py`. Stored data only (`results/fit_d357.json`, stored reference cells); no simulation. Output: `results/L_amend_gemini.json` (all numbers), `results/L_amend_gemini_preds.json` (per-replicate predictions). **Sanity assertion passed: re-derived registered R1-R5 hit counts equal the frozen `L_coverage.csv` in all 52 cells (0 mismatches).** **R4 stays as registered; its numbers below are unchanged from RESULTS.md.** R4', R2-lin and the Lambda projection are CORRECTIVE and post hoc in origin (suggested after the results were known); they cannot change the H1 verdict (R1 and R5). The stored truth cells are called reference estimates here.

## 1. Coverage and median prediction/reference (300 replicates; primary = absolute sigma; scaled = chi2/dof sensitivity)
Columns 1e-3 / 1e-4 / 1e-5 / 1e-6.

| recipe | d | coverage % primary | coverage % scaled | median pred/ref |
|---|---|---|---|---|
| R4 (registered, p^(d/2)) | 3 | 17.3 / 0 / 0 / 0 | 32.0 / 0 / 0 / 0 | 1.44 / 4.22 / 12.8 / 41.2 |
| R4 | 5 | 10.0 / 0 / 0 / 0 | 21.7 / 0 / 0 / 0 | 1.49 / 4.08 / 12.8 / 41.6 |
| R4 | 7 | 23.7 / 0 / 0 / 0 | 33.3 / 0 / 0 / 0 | 1.45 / 4.81 / 15.2 / 42.3 |
| **R4' CORRECTIVE (p^((d+1)/2))** | 3 | 95.0 / 95.7 / 95.3 / 95.0 | 96.0 / 96.3 / 96.3 / 96.0 | 0.99 / 1.02 / 1.00 / 1.01 |
| R4' | 5 | 95.0 / 95.0 / 95.0 / 94.7 | 95.0 / 95.3 / 95.3 / 95.0 | 1.03 / 0.99 / 0.99 / 1.02 |
| R4' | 7 | 94.0 / 82.7 / 81.3 / 93.0 | 95.0 / 87.0 / 86.3 / 94.3 | 1.00 / 1.17 / 1.18 / 1.04 |
| R2 (registered, constant prefactor) | 3 | 0 / 0 / 0 / 0 | 2.3 / 3.7 / 1.0 / 2.3 | 0.84 / 0.85 / 0.82 / 0.84 |
| R2 | 5 | 0 / 0 / 0 / 0 | 9.3 / 0.3 / 0.3 / 2.7 | 0.81 / 0.77 / 0.77 / 0.79 |
| R2 | 7 | 0 / 13.7 / 27.3 / 0 | 56.0 / 100 / 100 / 96.0 | 0.78 / 0.93 / 0.94 / 0.83 |
| **R2-lin CORRECTIVE, p^h(c0+c1 p)** | 3 | 91.3 / 93.7 / 91.3 / 92.3 | 94.0 / 96.0 / 94.0 / 95.0 | 0.98 / 1.01 / 0.98 / 1.00 |
| R2-lin | 5 | 68.3 / 84.0 / 83.7 / 65.3 | 73.0 / 88.3 / 87.7 / 69.3 | 1.07 / 1.05 / 1.05 / 1.08 |
| R2-lin | 7 | 6.0 / 0 / 0 / 0.3 | 16.3 / 0 / 0 / 1.7 | 1.18 / 1.45 / 1.47 / 1.30 |

(R4 medians reproduce RESULTS.md's 1.4/4.2/13/41-42x.) **The exponent error accounts for essentially all of R4's failure:** with the same fit and the same intervals, only the exponent changed and median prediction/reference falls from 41x to 1.01x at d=3 and d=5, and coverage goes from 0% to about 95% (nominal) at d=3 and d=5 at every target. At d=7 R4' is 93-94% at 1e-3 and 1e-6 but 81-83% at 1e-4/1e-5 (primary; 86-87% scaled) with median ratio 1.17-1.18. Caution for d=7: its reference has 7-8% se at 1e-4..1e-6, and R4' is flagged in section 4. R2-lin shows the opposite pattern: acceptable at d=3 (91-94%), degraded at d=5, near-zero at d=7 (median 1.3-1.5x high); adding a linear term to a fixed-exponent power law is not a general repair.

## 2. Is R4's 1.4/4.2/13/41x pattern explained by the p^(-1/2) exponent error? (registered rule)
- (a) median R4'/reference in [0.5, 2] in 12/12 (d,target) cells (range 0.99-1.18): **met**.
- (b) log share ln(median R4/R4') / ln(median R4/ref) at 1e-4, 1e-5, 1e-6, every d: d=3 0.98/1.00/1.00; d=5 1.01/1.00/0.99; d=7 0.90/0.94/0.99, all >= 0.80: **met**.
- **Verdict: EXPLAINED.** Mechanism (b, registered rule: all nine R4/R4' successive-decade ratios in [2.5, 4.0]): the ratios are 2.83, 3.13, 3.16 at every d (identical across d), against sqrt(10) = 3.16: **CONFIRMED**. The first step is below sqrt(10) (2.83) because the c1 p, c2 p^2 terms absorb part of the factor inside the window and the extrapolation is less than a full decade from the window centre. R4/ref successive steps are 2.7-3.3 (d=7: 3.32, 3.16, 2.78).

## 3. R5 decomposition (unrounded stored values)
Reference projection = p_L(5)^2/p_L(3) from the reference cells; distance factor = projection / d=7 reference; fit-input factor = R5 median / projection.

| target | ref d=3 | ref d=5 | ref d=7 | projection | distance factor | R5 median | R5/ref | fit-input factor |
|---|---|---|---|---|---|---|---|---|
| 1e-3 | 4.709e-4 | 6.108e-5 | 6.801e-6 | 7.924e-6 | 1.165 | 1.765e-5 | 2.59 | 2.23 |
| 1e-4 | 4.653e-6 | 6.432e-8 | 5.718e-10 | 8.892e-10 | 1.555 | 5.577e-9 | 9.75 | 6.27 |
| 1e-5 | 4.794e-8 | 6.435e-11 | 5.649e-14 | 8.639e-14 | 1.529 | 1.781e-12 | 31.5 | 20.6 |
| 1e-6 | 4.722e-10 | 6.244e-14 | 6.412e-18 | 8.256e-18 | **1.288** | 5.688e-16 | **88.7** | **68.9** |

ChatGPT's estimates verified: distance factor 1.17-1.56 (1.29 at 1e-6) and fit-input factor about 69 at 1e-6 (88.7 total). So about 94% of R5's log error at 1e-6 (ln 68.9 / ln 88.7 = 0.94; at 1e-4 0.81) comes from extrapolating the fitted rates in p, not from the distance step. R5 is a physical-error-rate extrapolation failure (the R1 exponent), not a distance-extrapolation failure.

## 4. Reference-uncertainty bounds (one shared reference value varying over ref ± 1.96 se)
Computed for all 152 cells (registered R1-R5 plus R4', R2-lin; primary and scaled). **34 cells can change verdict class** (miscalibrated <0.90 / calibrated 0.93-0.97 / other). Among the 52 registered primary cells, 47 cannot flip; the 5 flagged are all d=7 fixed-exponent: R2 d=7 at 1e-4 and 1e-5 (point 13.7%, 27.3%; hit fraction range 0-94%), R3 d=7 at 1e-4, 1e-5, 1e-6 (37.7%, 57.0%, 0%; ranges 0-95.7%, 0-95.7%, 0-92.7%). Reference relative se at those cells is 7-8% (time-capped runs), while the R2/R3 bias is 6-25%, so these five cannot be called miscalibrated or not from this data. **R1 and R5 cells (the H1 cells) and all R4 cells are not flippable: H1 is unaffected.** The "all 52 below 90%" statement therefore holds at the point estimates, and under reference uncertainty for 47/52 cells; 5/52 are undetermined. R4' d=7 at 1e-3, 1e-4, 1e-5 and 1e-6 (primary) is also flagged (e.g. 1e-4: 82.7% at the point, range 52.7-93.3%); R4' d=3 and d=5 cells are not flagged except R4' d=3 at 1e-5 under scaled sigma. Full per-cell table: `results/L_amend_gemini.json` key `reference_uncertainty`.

## 5. Single-operating-point Lambda projection (p = 4e-3, d=3 -> d=7, window point only)
P7hat = P5^2/P3 against the d=7 reference 1.709e-3: coverage **78.0%** (234/300; CP 72.9-82.6%), median prediction/reference 1.18 (IQR 1.07-1.27), median interval width ratio (upper/lower) 1.64. Biased high by about 18%, but this is a 4e-3 point, near the 8.9e-3 finite-size crossing's neighbourhood, and the single-point test is much less demanding than the far-extrapolation targets; it says the distance step alone is not exact even where data are plentiful. One cell, not a general statement.

## 6. d* table (ILLUSTRATIVE; p = 1e-3; target p_L = 1e-12; extrapolating in d beyond 7 is itself the assumption)
Rule as registered (Lambda_{5->7} from each recipe's own medians; R5: Lambda from R1 d=3,5). The circuits use rounds = d, so p_L is per d*-round memory experiment; the qubit-rounds column is the spacetime volume 2 d*^3.

| source | p_L(7) at 1e-3 | Lambda | d* | rounds | qubits 2d*^2 | qubit-rounds 2d*^3 |
|---|---|---|---|---|---|---|
| reference | 6.80e-6 | 8.98 | 23 | 23 | 1058 | 24,334 |
| R1 | 1.53e-5 | 6.64 | 25 | 25 | 1250 | 31,250 |
| R2 | 5.31e-6 | 9.28 | 21 | 21 | 882 | 18,522 |
| R3 | 5.44e-6 | 8.64 | 23 | 23 | 1058 | 24,334 |
| R4 (registered) | 9.85e-6 | 9.26 | 23 | 23 | 1058 | 24,334 |
| R4' | 6.77e-6 | 9.26 | 23 | 23 | 1058 | 24,334 |
| R2-lin | 7.99e-6 | 8.21 | 23 | 23 | 1058 | 24,334 |
| R5 | 1.76e-5 | 5.67 | 27 | 27 | 1458 | 39,366 |

Surprise for the reader: at 1e-3 the large factor errors barely move d* (R4 is 1.45x high yet picks the same d* = 23); errors in Lambda matter more (R1 +2 and R5 +4 distances, R2 -2). The 41x R4 miss is at 1e-6, a different design point; this table is at 1e-3 only.

## 7. Figures (`figs_amend/`; plotting only)
- `figA_prefactor.png/pdf`: p_L/p^((d+1)/2) vs p, 1e-6 to 1.5e-2; reference estimates, window medians, exact polynomials (d=3, d=5; surrogate with decoder frozen at 1e-3, d=7 incomplete and not drawn), R2's fitted prefactor (dotted; medians 394.9, 4.93e4, 5.31e6 for d=3,5,7). The prefactor drifts downward across the window (the fitted constant is a window average that sits 15-25% below the low-p plateau) which is exactly why R2 and R3 sit at 0.78-0.96 of the reference.
- `figB_money.png/pdf`: p_L vs p per d: window medians, reference estimates, R1 and R2 median 95% ribbons (the R2 ribbon is too thin to see, which is the over-confidence finding), finite-size crossing of the d=5 and d=7 window medians at p = 8.9e-3 (log-linear interpolation in the window; a marker, not a claimed threshold).
- Terminology: the reversal of the distance ordering near 8e-3 is a **finite-size crossing**.

The paper .tex was not edited.
