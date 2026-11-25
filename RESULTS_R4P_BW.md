# RESULTS: Amendment R4P-BW (R4' and R2-lin on the below-threshold window)
Pre-registered in PLAN.md (commit 18c2333) before any computation; CONFIRMATORY. Stored replicates `results/fit_bw.json` (300 reps, d=3,5,7, p=1.5-4e-3), same truth cells, targets, hit rule (primary absolute sigma), CP intervals, band [0.93,0.97]. Code `code/r4p_bw.py`; full output `results/L_r4p_bw.json`. Sanity: re-derived R1-R5 hits equal frozen `L_bw_coverage.csv` (0 mismatches). No simulation, no registered result changed.

## Bar: R4' d=3 and d=5 in the 93-97% band in >= 75% of cells (>= 6 of 8)
**FAIL: 5 of 8 (62.5%).** The three misses are all d=3 and all ABOVE 97% (over-covering): 1e-4 0.977, 1e-5 0.980, 1e-6 0.980. d=5 is 4/4 in band (0.943-0.963); d=3 at 1e-3 is in band (0.963). No cell below 90%. Median R4'/truth is 0.97-1.03 at d=3,5 at every target (R4 with d/2 was 2.4x, 7.3x, 23x at 1e-4, 1e-5, 1e-6): the exponent correction removes the bias; what fails the bar is slight over-coverage at d=3. (Sensitivity, chi2/dof-scaled: 4/8 in band, the other 4 above 97%: same direction, not better.)

## d=7 (no bar, separate)
Point-reference coverage 0.930, 0.947, 0.947, 0.950 (1e-3..1e-6): 4/4 in band, none below 90%, median ratio 0.94-1.09. **All four cells are UNDECIDED:** the verdict class flips within truth +/- 1.96 se (hit-fraction range 0.82-0.95, 0.84-0.96, 0.83-0.96, 0.92-0.96; d=7 reference se 3-8%). The d=7 reference is too uncertain to say R4' is calibrated or not there.

## Reported, not tested
| recipe (primary) | d=3,5 in band | d=3,5 below 90 / above 97 | d=7 in band | d=7 below 90 |
|---|---|---|---|---|
| R4' | 5/8 | 0 / 3 | 4/4 (all undecided) | 0 |
| R2-lin | 6/8 | 0 / 0 (2 cells at 0.927 fall in the 90-93 gap) | 1/4 | 3 (0.893, 0.707, 0.680) |
| R4 (d/2, registered) | 1/8 | 6 / 0 | 1/4 | 3 |
R2-lin scaled: 8/8 in band at d=3,5 (exploratory), but d=7 stays poor (1/4; 2 below 90). R2-lin d=7 at 1e-4 and 1e-5 is clearly miscalibrated at the point reference and over the whole reference range upper end (max 0.947).

## Reading
The registered R4 failure below threshold is the exponent: with (d+1)/2 the Bravyi-form recipe tracks the reference (median ratio ~1) and has no cell under 90% at any d. It misses the pre-registered bar only by over-covering at d=3 (3 cells at 97.7-98.0%); the bar as written is FAIL. Whether R4' is calibrated at d=7 cannot be decided with the current reference. R2-lin is calibrated at d=3,5 but not at d=7.
