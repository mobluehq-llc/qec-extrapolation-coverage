# RESULTS: Amendment 2, rolling-origin discrepancy inflation for R2/R3 (2026-11-24)
Pre-registered in PLAN.md Amendment 2 (commit 1b865c2) before any computation. [internal operations note redacted] Outputs: `results/L_inflate_coverage.csv`, `results/L_inflate_verdicts.txt` (all variants, all cells). Sanity: uninflated hits reproduce the stored coverage tables exactly (0 mismatches in both windows).

## Verdict (CONFIRMATORY, frozen rule): FAIL
Bar: coverage 93-97% in >= 18/24 R2/R3 cells in BOTH windows. Result: below-threshold (BW) window 5/24; straddling window 3/24.

| window | recipe | in band /12 | <90% | >97% | median cov (inflated; uninflated) | median width ratio (ln-scale; linear) |
|---|---|---|---|---|---|---|
| BW | R2 | 4 | 5 | 1 | 0.920; 0.657 | 2.21x; 2.22x |
| BW | R3 | 1 | 7 | 4 | 0.658; 0.063 | 2.84x; 2.85x |
| STRADDLE | R2 | 3 | 0 | 8 | 1.000; 0.000 | 14.7x; 16.0x |
| STRADDLE | R3 | 0 | 0 | 12 | 1.000; 0.000 | 21.8x; 23.3x |

- BW far targets (<= 1e-4): 3/18 R2+R3 cells in band, 11/18 still below 90% (R3 d=5 and d=7 remain the failures; median R3 coverage 0.66). Inflation raises coverage a lot (R2 median 0.66 to 0.92) but not into the band, and R3 stays miscalibrated.
- Straddling window: coverage becomes 100% in most cells by interval width (15-23x), not by calibration; widths exceed the triage's 10x useless-width criterion. The rule does not heal the unhealable in a calibrated way; it over-covers.
- Control R1 (side check, same inflation): BW 7/12 in band, median 0.942 (uninflated 5/12, 0.925), width 1.02x; straddling stays at 0.07 (not healed, as it should).
- EXPLORATORY (cannot carry the verdict): chi2/dof-scaled gives 4/24 (BW), 4/24 (straddle); leave-LOWEST-p-out gives 5/24 and 3/24 (BW R3 median cov 0.80, width 3.05x).

Registered caveat held: the one-step held-out residual underestimates discrepancy at far targets. No tuning was done after seeing results.
