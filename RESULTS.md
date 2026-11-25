# L RESULTS: extrapolation-interval coverage (2026-11-23; run exactly as frozen in PLAN.md FREEZE ADDENDUM)

[internal operations note redacted] **Reproducibility check:** all 4500 regenerated window points (k,n) are bit-identical to the 11-21 file (kept as `results/fit_d357_1121.json`; sha256 prefix 89f4f9322922c355 both). Analysis: [internal operations note redacted]

## Headline
**H1 SUPPORTED** (primary, and under the chi2/dof-scaled sensitivity): R1 is miscalibrated in 9/9 cells (CP-firm 9/9) and R5 in 3/3 cells (3/3 firm) at targets <= 1e-4. Actual coverage is not borderline: **0.0% in 36 of the 52 (recipe,d,target) cells, and all 52 are below 90%**, median coverage 0.000 for R1-R4 and 0.023 for R5. No recipe has a single cell in the 93-97% band (primary); under the scaled-sigma sensitivity R2 has one such cell (d=7, 1e-6, 96.0%) and R2 d=7 reaches 100% at 1e-4 and 1e-5 (intervals inflated by chi2/dof), the only form that ever over-covers.

Claims ledger: #1 (R1 <90% at 2-4 decades) 9/9 cells. #2 (R5 <90%) 3/3 (4/4 incl. 1e-3). #3 "at least one recipe within 93-97%, or none": NONE does (primary 0/52 cells). #4 truth validation passed 17/17 (11-21).

Direction (median prediction/truth over replicates; computed from the same fits): **R1, R4, R5 over-predict, increasingly with decades**: R1 d=7 2.2x at 1e-3, 8.6x at 1e-4, 28x at 1e-5, 80x at 1e-6 (d=5: 1.7/3.5/7.9/18; d=3: 1.2/1.9/2.8/4.3); R4 1.4x -> 42x (1e-3 -> 1e-6), same at every d; R5 2.6x/9.8x/32x/89x. **R2 and R3 (fixed exponent (d+1)/2) are close but slightly LOW: 0.73-0.96 of truth at every cell**, yet still 0% covered because their absolute-sigma intervals are only ~1-2% wide: a ~15-25% bias is 10+ sigma. So two different failures: free-exponent/Lambda-from-small-d forms are biased by orders of magnitude, fixed-exponent forms are nearly right but wildly over-confident.

## Coverage table (CP 95% in brackets), 300 replicates per cell
Truth = real-circuit IS (run1 + independent run2, inverse-variance); combined relative se: d=3,5 1.3-1.4%; d=7 3.2% (1e-3) and 7-8% (1e-4..1e-6, the 2400 s cap bound, 9-10% per run). Run1/run2 agree in all 12 cells (|z| <= 1.64). d=7 truth uncertainty (7%) is far smaller than the ~100x misses, but is not propagated, per the freeze.

**PRIMARY (absolute sigma)** coverage % (hits/300), columns = target p 1e-3 / 1e-4 / 1e-5 / 1e-6

| recipe | d | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|---|
| R1 | 3 | 60.0 [54.2,65.6] | 23.0 [18.4,28.2] | 18.3 [14.1,23.2] | 14.0 [10.3,18.4] |
| R1 | 5 | 2.0 [0.7,4.3] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R1 | 7 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R2 | 3 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R2 | 5 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R2 | 7 | 0.0 [0.0,1.2] | 13.7 [10.0,18.1] | 27.3 [22.4,32.8] | 0.0 [0.0,1.2] |
| R3 | 3 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R3 | 5 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R3 | 7 | 0.0 [0.0,1.2] | 37.7 [32.2,43.4] | 57.0 [51.2,62.7] | 0.0 [0.0,1.2] |
| R4 | 3 | 17.3 [13.2,22.1] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R4 | 5 | 10.0 [6.8,14.0] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R4 | 7 | 23.7 [19.0,28.9] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R5 | 7 | 5.7 [3.3,8.9] | 2.0 [0.7,4.3] | 2.0 [0.7,4.3] | 2.7 [1.2,5.2] |

**SENSITIVITY (chi2/dof-scaled)** coverage % (hits/300), columns = target p 1e-3 / 1e-4 / 1e-5 / 1e-6

| recipe | d | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|---|
| R1 | 3 | 68.0 [62.4,73.2] | 40.3 [34.7,46.1] | 35.3 [29.9,41.0] | 28.7 [23.6,34.1] |
| R1 | 5 | 16.0 [12.0,20.6] | 5.7 [3.3,8.9] | 3.0 [1.4,5.6] | 2.0 [0.7,4.3] |
| R1 | 7 | 10.3 [7.1,14.3] | 1.7 [0.5,3.8] | 1.3 [0.4,3.4] | 1.7 [0.5,3.8] |
| R2 | 3 | 2.3 [0.9,4.7] | 3.7 [1.8,6.5] | 1.0 [0.2,2.9] | 2.3 [0.9,4.7] |
| R2 | 5 | 9.3 [6.3,13.2] | 0.3 [0.0,1.8] | 0.3 [0.0,1.8] | 2.7 [1.2,5.2] |
| R2 | 7 | 56.0 [50.2,61.7] | 100.0 [98.8,100.0] | 100.0 [98.8,100.0] | 96.0 [93.1,97.9] |
| R3 | 3 | 78.3 [73.2,82.9] | 87.7 [83.4,91.2] | 52.0 [46.2,57.8] | 75.0 [69.7,79.8] |
| R3 | 5 | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R3 | 7 | 1.7 [0.5,3.8] | 100.0 [98.8,100.0] | 100.0 [98.8,100.0] | 50.0 [44.2,55.8] |
| R4 | 3 | 32.0 [26.8,37.6] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R4 | 5 | 21.7 [17.1,26.8] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R4 | 7 | 33.3 [28.0,39.0] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] | 0.0 [0.0,1.2] |
| R5 | 7 | 28.3 [23.3,33.8] | 13.7 [10.0,18.1] | 14.7 [10.9,19.2] | 16.7 [12.6,21.4] |

## Controls (run as registered)
- **Control A (data-only noise, `results/L_control_a.txt`):** the data-only dataset is rejected by every recipe in 15/15 cells (reject rate 1.000 for R1-R4), i.e. it is detected as a different model. Fresh circuit data are rejected at 0.185 (R1), 0.490 (R2), 0.520 (R3), 0.062 (R4); no cell >= 0.95 for R1/R4, 2/15 for R2 and R3. So the fit-vs-fresh-data test itself has a false-reject rate well above 5% for R1-R3 on absolute sigma (the same overdispersion the scaled sensitivity addresses). Reported, not tuned.
- **Control B (edge-iid surrogate, `results/L_control_b.txt`):** all five recipes are miscalibrated in every cell (R1-R4 12/12, R5 4/4; median coverage 0.000-0.062), as the F2 positive control is registered to behave.
  Caution for reading the main table: B shows that even a model whose structure is simpler than the circuit gives ~0 coverage under these recipes, so the main result is not unique to circuit-level detail; the failure is in the extrapolation forms plus narrow absolute-sigma intervals.

## d=9 holdout (ONE SHOT, run once, registered in PLAN.md before the run)
Window data `fit_d9.json` (300 reps) and six independent truth runs (`truth.py 9 1e-4`, seeds 9001..9006, 3000 s cap each, 14.2-14.6 M samples each), scored by `code/holdout.py` unchanged.
- **Mechanical registered verdict (script output): "CONFIRMED (both <90%)":** R1 coverage 0/300 CP[0,1.2%], R5 0/300 CP[0,1.2%]; also R2, R3, R4 0/300. Median prediction/point-truth: R1 184x, R2 9.8x, R3 10x, R4 45x, R5 261x [CORRECTED 2026-11-24: these ratios are OVERSTATED and must not be quoted; the point truth 5.5e-13 lies below the certified lower bound 6.5e-13 (A-0). Against the certified interval the median R1 prediction is 6.3x (vs U) to 155x (vs L) high, R5 9.0x to 221x; see the CERTIFIED p_L INTERVAL section below and paper v5] (all over-predict the point estimate; note R2/R3 were 0.73-0.96x at d<=7, so their d=9 ~10x may reflect truth uncertainty or a real break, undetermined).
- **But the truth is NOT resolved, and this must travel with the verdict.** Five of six runs put p_L(d=9, 1e-4) at 3.4-4.8e-13 (mutually consistent); the sixth hit a rare heavy sample and gives 2.7e-11. The estimator's conservative se (zero-failure/low-ESS inflation on the high-weight strata w=9..12 with only 5000 samples each) is 3.4-3.8e-11 per run (one run 2.2e-10), i.e. relative se about 70-110x the estimate; the inverse-variance combination is 5.5e-13 +- 1.6e-11 (rel se 2905%). **That fails the method's own width bar (<= 20%) by two orders of magnitude.** The mechanical hit test used the point value; the coverage 0 is therefore informative only in direction. An interval-aware read: the conservative upper bound of the truth (~3e-11 at +1.96 se) is still below the median R1/R5 predictions (~1e-10), so the direction (over-prediction) survives, but the magnitude (10x-260x) does not.
- Registered handling: runs hit the cap, so per the freeze they are reported as-is and NOT extended or rerun. **Reading for the paper:** holdout direction consistent with H1; magnitude not established; d=9 at 1e-4 (p_L ~ 1e-12 or lower) is beyond what this IS method resolves in 50 min x 6 cores. A tighter truth needs a new, labelled amendment (more weight-9..12 samples), not a rerun under this label.

## Not done / limits
- One window, one noise model, one decoder-weight convention (same-p DEM). Per-d/per-window breakdowns beyond the tables are exploratory and were not run.
- Vuong comparison (exploratory) not run in this pass.
- Truth se for d=7 not propagated into the hit test (frozen).
- numpy portability: INTERESTING row 36 says seeded results differ across numpy versions; [internal operations note redacted]

## Below-threshold amendment (2026-11-24; CONFIRMATORY; [internal operations note redacted]
[internal operations note redacted] frozen code hashes matched; `results/fit_bw.json` sha256 prefix 4db811282c8030a8). Same R1-R5, same 12 real-circuit truth cells, same bars. Analysis `code/analyze_bw.py` (verdicts in `results/L_bw_verdicts.txt`, tables `results/L_bw_coverage*.csv`); the BW-2 ratio and the cross-table counts come from `code/bw_extra.py` and a csv comparison with `results/L_coverage.csv` (read-only, no new fitting). One slip: [internal operations note redacted] no coverage had been computed; fixed by copying the frozen truth files.

**Headline: the registered prediction BW-1 FAILED. Below threshold, H1 is NOT SUPPORTED for the free-exponent form.** Per the amendment's own rule this means: for R1 the H1 failure was a straddling-window artefact. R5 stays miscalibrated.

| prediction (registered) | result | verdict |
|---|---|---|
| BW-1: R1 and R5 each miscalibrated (<90%) in >= 2/3 of far cells | R1 0/9 (CP-firm 0/9), R5 2/3 (CP-firm 1/3); chi2-scaled: R1 0/9, R5 1/3 | **FAILED (R1)**, met for R5 only under primary sigma |
| BW-2: R1 coverage strictly higher than straddling table in >= 6/12 cells; d=7, 1e-6 over-prediction closer to 1 than 80x | 12/12 higher; d=7 1e-6 median prediction/truth 0.85 | PASSED |
| BW-3: <= 2 of 16 R1/R5 cells in the 93-97% band; R1 median coverage < 0.90 | 5/16 in band (all R1); R1 median 0.925 | FAILED |
| BW-4: R2 and R3 miscalibrated in >= 2/3 of cells | 19/24 cells (R2 11/12, R3 8/12); far set 16/18 | PASSED |

Primary (absolute sigma) coverage, % (hits/300), columns 1e-3 / 1e-4 / 1e-5 / 1e-6:
R1 d3 93.0/91.0/91.3/91.3; d5 92.0/93.7/94.0/92.7; d7 92.3/92.3/94.0/94.0 (median 92.5; all 12 cells between 91.0 and 94.0).
R2 d3 68.0/83.7/35.0/63.3; d5 94.7/33.3/32.3/79.3; d7 85.7/0.0/0.0/72.0.
R3 d3 96.0/95.0/84.3/94.0; d5 0.0 in all four; d7 94.7/0.0/0.0/12.7.
R4 d3 96.0/2.7/0.0/0.0; d5 91.7/6.0/0/0; d7 94.0/3.0/0/0.
R5 d7 71.3/82.7/88.3/91.3.
Chi2-scaled sensitivity: R1 10/12 cells in 93-97% (median 94.3%), R5 median 88.7%.

Median prediction/truth (1e-3, 1e-4, 1e-5, 1e-6): R1 d3 1.00/1.15/1.27/1.45, d5 1.03/1.11/1.26/1.47, d7 0.95/1.08/1.02/0.85 (straddling window: 2.2/8.6/28/80x at d=7). R2 0.94-1.20, R3 0.88-1.21 (no longer systematically low: the straddling-window "close but low" pattern did not persist). R4 over-predicts 2.3-2.6x at 1e-4 up to 22-24x at 1e-6 at every d (unchanged failure of the quadratic-in-p form). R5 1.26-1.92x.

**Reading.** (1) On a window wholly below the crossing, the plain free-exponent log-log fit with its absolute-sigma 95% interval is close to calibrated at 1-3 decades (91-94% in 12/12 cells), so the main-table headline "R1 fails by orders of magnitude" is a property of fitting through the threshold, not of free-exponent fits in general. (2) What does not go away: R4 (24x), R5 (Lambda from d=3,5: 1.3-1.9x over, 71-91% coverage, falling at the 1e-3 end), and R2/R3's over-confident fixed-exponent intervals (R3 d=5 0% in every cell; R2 d=7 0% at 1e-4, 1e-5). Fixed-exponent forms are now biased by only 0-20% but their intervals are about 1-2% wide. (3) Limits: one code, one noise model, one decoder-weight convention; the 1e-3 target is only 0.18 decades from the window, so the far cells here are 0.2-2.2 decades below the window versus 0.6-2.6 before (stated in the amendment); the d=7 truth se (7-8% at far cells) is not propagated and is comparable to the 91-94% R1 coverage margin, so "calibrated" here means "not rejected at this truth precision". The original table, controls and d=9 holdout are unchanged and still describe the registered straddling window.

### Notes added with paper v4 (2026-11-24; no result changed)
- **d=7 runs and the time limit.** Checked the stored run records (`results/truth_d7_p*.json`, `truthextra_d7_p*.json`: [internal operations note redacted] The estimator ran 12 times at d=7: 9 ended at the 2400 s limit (all 8 runs behind the four reference cells, plus the validation run at p=2e-3, rel se 2.1%), 1 ended at the 4e7-draw limit after 2148 s (p=1.5e-3, rel se 2.4%), 2 stopped at the 2% rule (p=3e-3, 4e-3). "Every run at d=7 reached its time limit" (paper v3) is therefore true of the reference runs only; v4 says so.
- **Decades below the window.** The amendment text and the Reading above give the far set as 0.18..2.2 (new) and 0.6..2.6 (old) decades below the lower window edge. Those are one decade too small: for targets 1e-4..1e-6 the correct figures are 1.2..3.2 and 1.6..3.6. The difference between windows (0.43 decades) is unchanged; no rule or result uses these figures.
- **Interval widths below threshold** (`paper/make_numbers_v4.py`, descriptive, same fits): cell-median 95% half-width of R1 is 15% / 61% / 126% / 219% of the prediction at 1e-3 / 1e-4 / 1e-5 / 1e-6; R2 5.1%, R3 2.9-4.6% (not "1-2%", which is nearer their standard error). Median fitted R1 exponents 1.95 / 2.95 / 4.03 (straddling window 1.82 / 2.65 / 3.49). 4 of 12 R1 cells have a CP upper limit below 95%.
- **Not run:** R4' (corrective), controls and the reference-uncertainty bounds on the below-threshold window; they were not in the amendment.

---
# CERTIFIED p_L INTERVAL (2026-11-24; [internal operations note redacted] pre-registered in PLAN.md "Certified interval amendment")
Method: exact stratum masses P(W=w) x per-stratum betting confidence bounds on the IS mean (rigorous range bound b_w = B_w/0.3), delta=0.05 union-bounded. Assumptions A1-A5 are listed in PLAN.md; A-0 (f_w=0 for w<=4 at d=9) is an ASSUMPTION, supported but not proven.

**Validation (d=3/5, exact truth, 24 replicates per cell, run as registered): ALL PASS.** Truth inside [L,U] 24/24 in every cell. Median relative width: d=3 p=1e-3 2.0%, 1e-4 1.5%; d=5 p=1e-3 14.1%, 1e-4 11.0% (bar 50%). Data: `results/cert_validation/`.

**d=9, p=1e-4 (NEW measurement, not the holdout rerun; 1.0e9 fresh draws over w=5..12 + 1e7 each for w=1..4; [internal operations note redacted] chunks in `results/cert_d9_chunks/`, result `results/cert_d9_result.json`):**
- **Under A-0: p_L in [6.5e-13, 1.6e-11] (95%, certified), U/L = 25; plug-in point estimate 3.8e-12.** Dominant terms: w=5 (P=1.0e-4, 37,126 failing draws in 8.8e8) contributes [5.4e-13, 8.6e-12]; w=6 [8.7e-14, 5.1e-12]; w=7 U 1.7e-12; tail P(W>12)=3.6e-15.
- **Assumption-free (strata 1-4 bounded statistically): [6.3e-13, 1.1e-6]: VACUOUS**, as pre-registered. Zero failures in 1e7 draws cannot certify below ~3e-6 per stratum; w=1,2 alone contribute 8.7e-7 and 2.0e-7 to U. Everything below 1e-6 therefore rests on A-0.
- Consequence for the old holdout: its point truth 5.5e-13 (rel se 2905%) lies just BELOW the certified lower bound 6.5e-13: the holdout under-estimated (missed the rare heavy-weight failing sets, consistent with its 6th run at 2.7e-11). The magnitudes "R1 184x, R5 261x" in the holdout section are therefore overstated; against the certified interval R1's median prediction (1.0e-10) is 6.3x above U and 155x above L (paper v5, make_numbers_v5.py keys c9.pred.overU.R1 / c9.pred.overL.R1; the earlier "about 6-100x" here was a rough figure), not 184x.
- **Decisive-for-R1/R5 bar (>= 285/300 replicates, both, under A-0): NOT MET. Not decided jointly.** Certified missed (U < lower end of the recipe's 95% interval): **R1 300/300** (median interval [6.4e-11, 1.6e-10]), **R4 239/300**, **R5 208/300** (R5 intervals are wide, median [2.6e-11, 7.9e-10], so U=1.6e-11 sits below the lower end in only 69% of replicates; 0/300 above), R2 0/300 and R3 0/300 (their intervals, ~5.1-5.8e-12, lie INSIDE [L,U]: not rejected, not confirmed). Under assumption-free: 0/300 for all recipes (vacuous).
- So: given A-0, the certified d=9 interval is tight enough to say **R1 missed (over-predicted, 300/300)**; R5 missed in 69% of replicates, short of the registered 95% bar; R2/R3 (fixed-exponent forms) are consistent with the certified interval. The H1 holdout direction survives for R1; for R5 it is only partly certified.
- Not done: d=7 interval (not needed to answer the question; weights 1-3 exact, 0 failing sets, would support A-0 but was not run); anytime-valid stopping, R6 and other triage candidates; prior-art search (triage: [internal operations note redacted] the empirical bar for A-0 at d=9 (only 1e7 draws per w<=4 with 0 failures plus d=7 exhaustive evidence).
