# Checks of the Gemini and ChatGPT red-team claims on L v5 (2026-11-24)

Rules were fixed before computing in `PLAN.md` ("CORRECTIVE AMENDMENT (Gemini red-team v5 checks)" and its addendum; ERRATA.md entry written first). Stored data only; no simulation. [internal operations note redacted] Sanity assertions passed: re-derived R5 and R2 d=7 hit counts equal the frozen `L_coverage.csv` and `L_bw_coverage.csv` (0 mismatches, both windows); IVW reference reproduces the frozen truth to 1e-12; frozen block-rate d* reproduced for all 8 sources. The .tex is not edited. Everything here is CORRECTIVE (post hoc in origin) and changes no frozen verdict.

## 1. Per-round normalisation: Gemini is right in substance, and more changes than the 25/21
Premise checked in code: `rounds=d` (`truth.py`, `fitdata.py`), so p_L(d) is a block rate per d-round experiment, and R5 and Table VII apply a geometric step to it. Identity: if eps(d) = (1-(1-2p_L)^(1/d))/2 is geometric, p_L(5)^2/p_L(3) = (25/21) p_L(7) exactly.

Reference-value distance factor (Table VI, "Distance"), block versus per-round:

| p | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|
| block (paper) | 1.165 | 1.555 | 1.529 | 1.288 |
| per-round | **0.978** | **1.306** | **1.285** | **1.082** |
| Lambda per-round 3->5 / 5->7 | 12.85 / 12.57 | 120.6 / 157.5 | 1242 / 1595 | 12605 / 13634 |

So at 1e-3 the per-round suppression is constant to 2% (0.98) and at 1e-6 to 8%; 1.3 remains at 1e-4 and 1e-5, which is where the d=7 reference sits lowest (item 2). The claim "the whole 1.17-1.29 is the normalisation" is right at 1e-3 and 1e-6, not at 1e-4/1e-5. Gemini's statement that Table VI's entire range is the bug is too strong; "about 1.19 of it" is exact.

R5 with the per-round projection (each interval shifted by the replicate's own factor, which is 0.840 = 21/25 to 3 digits at these p_L):

| window | p | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|---|
| first (hits /300) | block | 17 | 6 | 6 | 8 |
| | per-round | 44 | 8 | 8 | 8 |
| | median R5/ref block -> per-round | 2.59->2.18 | 9.75->8.19 | 31.5->26.5 | 88.7->74.5 |
| second (coverage) | block | 71.3% | 82.7% | 88.3% | 91.3% |
| | per-round | 93.7% | 88.7% | 90.3% | 93.3% |
| | median R5/ref | 1.26->1.06 | 1.72->1.44 | 1.92->1.62 | 1.89->1.59 |

H1 (first window, far cells 2.0-2.7%) is unaffected. On the second window R5's registered "miscalibrated at >= 2/3 of far cells" (1e-4 and 1e-5 were below 90%) would become 1 of 3 (88.7% only; 90.3% is above the line by 0.3 points); BW-1 had already failed through R1, so its verdict does not change, but the sentence "R5 is mixed ... 71% rising to 91%" would read 94% / 89% / 90% / 93% under per-round.

d* table (p = 1e-3, 1e-12 per d*-round memory experiment, per-round Lambda from the same medians), block -> per-round:

| source | d* block | d* per-round | qubit-rounds per-round |
|---|---|---|---|
| reference | 23 | **21** | 18,522 |
| R1 | 25 | 23 | 24,334 |
| R2 | 21 | 21 | 18,522 |
| R3 | 23 | 21 | 18,522 |
| R4 | 23 | 21 | 18,522 |
| R4' | 23 | 21 | 18,522 |
| R2-lin | 23 | 21 | 18,522 |
| R5 | 27 | 23 | 24,334 |

Consequences for the text: "R2 asks for d*=21, two fewer than the reference, which would under-provision" no longer holds (R2 equals the reference); R5's inflation is 24,334 versus 18,522 (1.31x), not 39,366 versus 24,334 (1.62x); the statement that R4 hides a 1.45x miss inside the same d* as the reference remains true (all the fixed-exponent fits match). ChatGPT's remark that Table VII's arithmetic works is correct for the block-rate rule as stated; the rule is what is questionable. Recommendation: present Lambda per round, keep the block-rate numbers in the supplement as the registered version. ERRATA.md has the entry.

## 2. The d=7 reference "dip": mostly not a weighting artefact, and the premise is wrong
(a) Combination. Both runs per cell have almost the same se, so inverse-variance weights are 0.50/0.50 at 1e-3 and 1e-4, 0.48/0.52 at 1e-5, 0.61/0.39 at 1e-6. Equal weights change the combined value by 0.00%, 0.00%, +0.33%, +0.80% at 1e-3..1e-6. Run-to-run differences are z = 0.41, 1.19, 1.02, -0.49. Pooling is not what produces the 15% dip; it is in both runs at 1e-4 and 1e-5 (run 1: 6.19, 6.12; run 2: 5.24, 5.22 in units of 1e6 p^4 at 1e-4, 1e-5; at 1e-6 6.22 and 6.70; at 1e-3 6.89 and 6.71). Whether each run is individually biased low by missed heavy sets cannot be told from two seeds that share the estimator; the d=9 case (certified lower bound above the pooled point) shows the hazard is real for that estimator in that regime.
(b) Shape. Gemini's premise that p_L/p^4 "must approach c_4 monotonically" is not right as stated: at leading order p_L/p^4 includes the no-further-fault factor exp(-lambda(p)) (lambda = sum of mechanism probabilities, 0.2106 at 1e-4, 0.0211 at 1e-5, 0.0021 at 1e-6), which makes the ratio rise as p falls until about 1e-5, plus positive higher strata that make it larger at 1e-3. Expected 1e-4/1e-6 ratio from the factor alone: 0.812; observed 0.892 +- 0.091 (z = +0.88). Expected 1e-5/1e-6: 0.981; observed 0.881 +- 0.095 (z = -1.05). The pre-registered smooth model A exp(-lambda)(1+Bp) on all 8 runs FAILED (chi2 44.3 on 6 dof; the 1e-3 points need more than one linear term) and is not used. Post-hoc fallback, labelled so: A exp(-lambda) fitted to the six runs at p <= 1e-4 gives A = 6.38e6, chi2 = 6.4 on 5 dof (p = 0.27), largest |z| 1.8. So the 1e-4..1e-6 values are consistent with a smooth curve at their stated se; a downward bias of the size claimed is neither shown nor excluded (it would have to be about 10% at 1e-4 and 1e-5 and the se is 7-8%).
(c) R2 d=7 coverage re-scored (hits/300; 1e-3, 1e-4, 1e-5, 1e-6). First window: IVW (frozen) 0, 41, 82, 0; equal weights 0, 41, 69, 0; run 1 only 0, 0, 0, 0; run 2 only 0, 275, 266, 0; smoothed (post hoc) 0, 249, 0, 0. Second window: IVW 257, 0, 0, 216; equal 257, 0, 0, 247; run 1 202, 61, 25, 88; run 2 284, 0, 0, 284; smoothed 257, 0, 104, 181. So the R2 d=7 verdicts at 1e-4 and 1e-5 move between 0% and 92% depending on which run is the reference: they are undecided, as the paper already says (five registered cells), and the equal-weight change is small (at most 31 hits of 300, at 1e-6 on the second window; 13 on the first). Gemini's specific claim that the second-window 0.0% cells at 1e-4 and 1e-5 are an artefact of the dip is partly supported at 1e-5 (smoothed reference: 35%) and not at 1e-4 (0% in 4 of 5 variants; R2's median prediction is 1.16-1.29 of any of these references). Direction of the predictions: R2 is high there under every variant.
Effect on item 1: with the smoothed d=7 references the per-round distance factor would be 1.45 at 1e-4 and 1.16 at 1e-5 (not flat), i.e. the factor depends on the reference at the level of 15%.

## 3. d=7 reference refinement: pre-registered, not started
Rule (PLAN amendment item 3): fixed allocation from the existing 8 runs as a discarded pilot, x1.5 safety, floor 1e5 draws, per-stratum sums pooled, plain se for the 1.5% bar plus the conservative se reported, no 1/sigma^2 weights. [internal operations note redacted] draws needed for 1.5% plain se: 2.1e8 (1e-3), 1.13e9 (1e-4), 1.36e9 (1e-5), 1.31e9 (1e-6); strata 4 and 5 take more than 99% of the allocation at p <= 1e-4. Total 4.0e9 draws, 6.0e9 with the safety factor. Measured rate 1.0-1.4e4 draws/s/core (median 1.1e4): **about 100 CPU-h, 150 CPU-h with safety** (Gemini's "about 7 CPU-h per cell" is 4-5 times too low, because its 10x-the-budget rule assumes se falls as 1/10; the 2400 s runs gave 7-11% and se goes as 1/sqrt(draws): 7% -> 1.5% needs 22x). [internal operations note redacted] Entered in MACHINE_PLAN.md. A note on design: f_4 at d=7 is nearly p-independent (same decoder weights), so a shared stratum-4 pool across the three low p could cut the cost by about 3x; not adopted because it ties the four reference cells together.

## 4. SM Table S16 arithmetic: consistent
Sums reproduced from `cert_d9_result.json`: sum of contributions to L = 6.522e-13; sum to U with the tail P(W>12) = 3.63e-15: 1.6069e-11 (matches U; Gemini's 1.6064e-11 omits the tail and rounding). For the betting upper limit with n = 2000, S = 8 strata (so delta' = 0.05/16 = 0.003125): the failing draws carry importance scores of at most 9.1e-9 (w=11, 5 draws, mean 5.0e-9) and 4.0e-7 (w=12, 10 draws, mean 7.8e-8) against a range bound b = 3.40, so they do not move the bound: the limit computed with all 2000 draws set to zero is 0.0146335000277 (w=11) against the stored 0.0146335000401, and 0.014688783940 against 0.014688784332 (w=12): ratios 1.0000000008 and 1.00000003. Five versus ten failing draws are therefore irrelevant, and f11+ and f12+ are nearly equal because they are the zero-data bound for n = 2000 and b = 3.40 versus 3.41. Gemini's scale estimate b ln(1/delta')/n = 0.0098 is the single-lambda limit; the 18-value grid costs a further ln 18 = 2.89 in the exponent: b (ln 320 + ln 18)/n = 0.0147 (w=11) and 0.0148 (w=12), matching 0.01463 and 0.01469. Same for w=10 (0.001572 against 0.001583 predicted). So the bound is valid and conservative here; the point to note in the text is that for w >= 11 it is the zero-failure bound, and the 5 and 10 failing draws are uninformative because their weights are tiny. The 18 grid values: item (c).

## 5. Crossing wording
From `fit_d357.json` by the frozen log-linear rule: window means cross at 7.74 (d=3/5), 8.28 (d=3/7), 8.89 (d=5/7) x 1e-3; medians at 7.74, 8.30, 8.93. The text uses 8.28 (Sec. II.B, the d=3/7 crossing of the window means that sets the second window) and "about 8.9" (Fig. 1 and Sec. III.A, d=5/7 medians). The three-value statement Gemini attributes to SM S5 is not in the v5 supplement source (the three values are in PLAN.md only). Suggested wording: "the pairwise crossings of the window means are at 7.74, 8.28 and 8.89 x 1e-3 (d = 3/5, 3/7, 5/7); the window rule uses the d=3/7 value".

## (a) Sliding-window sweep (ChatGPT): pre-registered, not run
See PLAN addendum (a). Cost from stored shot counts and measured rates: d=3 0.02-0.06 CPU-h, d=5 5.7-17 CPU-h (windows), references up to 17 CPU-h, total 6-35 CPU-h; [internal operations note redacted] Design check: the cost model reproduces the stored second window at d=5 (0.35 CPU-h per 100 replicates against 0.23 measured).

## (b) "Exact" surrogate checks against a truncation: ChatGPT is right, and only partly certifiable
The surrogate's "exact" value is the polynomial through order 4. The manifests also hold, per complete shell w <= 4, the sums T_j for j <= 3 (orders up to w+3), and Poisson-binomial P(W >= 5) is exact. That gives a certified bracket for the full surrogate p_L (rule in PLAN addendum (b)):

| cell | 4th-order poly | certified [L, U] | known orders 5..7 of shells <= 4 relative to poly | bracket width / poly |
|---|---|---|---|---|
| d=3, 1e-3 | 4.6588e-4 | [4.6575e-4, 4.6620e-4] | -0.026% | 0.095% |
| d=3, 1e-4 | 4.7376e-6 | [4.737588e-6, 4.737592e-6] | -3e-5% | 1e-4% |
| d=5, 1e-4 | 6.2928e-8 | [6.2768e-8, 8.0179e-8] | -0.25% | 27.7% |
| d=5, 1e-3 | 6.2311e-5 | [5.19e-5, 1.0609e-3] | **-16.7%** | **16 (vacuous)** |

Re-scoring: the 96 certified replicate intervals contain the 4th-order polynomial in 96/96; the whole certified bracket lies inside [L,U] in 48/48 at d=3 (both p), and in 0/48 at d=5 (the bracket there is wider than the intervals; every interval does intersect it: 48/48). [internal operations note redacted] inside the 4th-order polynomial 19, 20, 20, 19 (d=3 1e-3, 1e-4; d=5 1e-3, 1e-4) = 78/80, as the paper says; all 78 also touch the certified bracket, so at d=3 the validation is against the full value, and at d=5 it rests on the polynomial plus the bracket's weak statement. The result worth reporting: at d=5, p=1e-3 the known higher-order terms of shells <= 4 alone lower the polynomial by 16.7%, so the "exact" value is not exact to better than the unknown shell >= 5 contribution (positive, at most 1.0e-3); the IS estimates agree with the polynomial to 1-2%, which means the unknown positive part nearly cancels the known negative part there. That is consistent but it is a coincidence the text should not call exact. Recommended wording: "exact through fourth order; the full value is certified to 0.1% at d=3 and not certified at d=5, p=1e-3". Closing it needs the shell-5 sum: d=3 is C(286,5) = 1.5e10 decodes (a few CPU-h, feasible); d=5 is 2.4e13 (not feasible).

## (c) The tail W > 120 and the grids
Chernoff bound for the Poisson-binomial W of mean lambda: P(W >= 121) <= exp(-lambda)(e lambda/121)^121 = 10^-241.1 at d=9, p=1e-4 (lambda 0.4545) and 10^-161.3 at d=7, p=1e-3 (lambda 2.104). Immaterial; add it to the sentence "strata beyond 120 are neglected". Betting grids (lambda_j = g_j / b_w, b_w the stratum's range bound; the product is averaged over the 18 values):
- Lower limit ("mean > m" test, prod(1+lambda(Y-m))): g = 0.001, 0.00183, 0.00336, 0.00617, 0.01131, 0.02074, 0.03803, 0.06975, 0.1279, 0.2346, 0.4301, 0.7888, 1.4465, 2.6527, 4.8645, 8.9207, 16.359, 30. A value is used only where lambda m < 1, i.e. m < b_w/g: all m for g <= 0.7888, m < 0.69 b_w, 0.38 b_w, 0.21 b_w, 0.11 b_w, 0.061 b_w, 0.033 b_w for the last six.
- Upper limit ("mean < m" test, prod(1+lambda(m-Y))): g = 0.001, 0.0015, 0.00225, 0.00338, 0.00508, 0.00762, 0.01145, 0.01718, 0.0258, 0.03873, 0.05814, 0.08728, 0.1310, 0.1967, 0.2953, 0.4433, 0.6655, 0.999; admissible for every m in [0, b_w) since g <= 0.999 gives lambda (b_w - m) < 1.
A bound is stated when the average reaches 2S/delta (=1/delta'); grids are fixed in advance, not tuned.

## (d) "registered" in S13
In the supplement source the caption of Table S13 (EXPLORATORY, reference se added in quadrature) ends "The bounds of Table S12 are the registered treatment of the same question", and Table S12 (the shared-reference bounds) is captioned CORRECTIVE and sits in a section titled "Corrective analyses (registered window)", where "registered" modifies "window" (the first window), and inside the caption it refers to the registered recipes R1-R5. Fix: S13 should read "the corrective treatment"; and ChatGPT's chronology point stands as worded in the paper's own freeze paragraph ("frozen before coverage was computed, after initial data generation"), which already says that. The 'S13' section Gemini and ChatGPT cite does not exist as a section number in the v5 source (12 supplementary sections); it is the table number.

## (e) Literature statements (labels now in PRIORART.md)
Mayer et al. exists (arXiv 2509.13678; first submitted Sep 2025, revised Oct 2026; splitting at circuit level, rates below 1e-20 claimed). Mullan, Weippert and Brown exists (arXiv 2607.27153, 29 Jul 2026) and has the d=11 hook example: a fit to the five highest-p points gives 9.47e-9 at p = 3e-4 against 1.54e-8 measured. Its pruning shrinks known failing patterns (worked example at d=7); it cannot prove absence. Beverland et al. A.3 / Fig. 16 fit p^(d/2) exp(a + bp + cp^2) directly to Monte Carlo for the Gross code with the Relay decoder (printed parameters [8, 11.43, 1881, -10303]) and compare with their failure-spectrum ansatz: same algebraic form as R4. All three read via a summarising fetch tool, so quotations are second-hand; the full papers should be read before citing.

## (f) A-0 on the real decoder: plan only
PLAN addendum (f): [internal operations note redacted] the pruning of Mullan et al. is a counterexample search, not a proof. Cost cannot be pinned before the pilot (pilot <= 20 CPU-h; the weight-3 enumeration fallback was 363-727 core-h and has no route for weight 4). Not started.
