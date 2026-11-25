# How often do extrapolated logical error rates contain the truth?

Reproducibility bundle for paper L (arXiv version v7): a coverage benchmark of 95% intervals from extrapolated logical error rates, for a rotated surface-code memory at d=3, 5, 7 under circuit-level depolarising noise with matching decoding. **Staged 24 November 2026; not yet published (no repository URL or DOI yet).**

**Author:** Adam Faust, MOBLUEHQ, Inc. See `CITATION.cff`.

## What the paper is
Resource estimates for fault-tolerant machines rest on logical error rates too small to simulate; they come from fits at high physical error rate p. The paper fixes the truth with reference estimates (fault-stratified importance sampling, exact low-weight enumeration, and a certified d=9 interval from betting-type bounds), fits five registered recipes R1-R5 (free power law, fixed exponent (d+1)/2 per distance, shared curve, smooth prefactor, extrapolation in distance) to 300 replicate data sets on two windows of p, and counts how often the 95% interval contains the reference. Registered, corrective and exploratory results are labelled as in the paper.

## One command
```
make verify          # python3 -I verify.py CLAIMS.md ; standard library only, about a minute
make verify-full     # additionally regenerates every number and table; needs numpy + scipy (pip install -r requirements.txt)
```
`make verify` (1) **re-fits R1-R5 on the stored replicates** (`results/fit_d357.json`, `results/fit_bw.json`: 300 replicates x d=3,5,7 x five window points) in pure Python, requires the hit counts of all four coverage tables (primary and chi2/dof-scaled, both windows; 52 cells each) to equal the stored CSVs, and recomputes the headline counts (24 of 24 and 19 of 24 fixed-exponent cells below 90% coverage; R2 and R3 upper-limit excesses; R1 coverage on the second window) and compares them with the paper's macros; (2) checks `CLAIMS.md`: every abstract number against its macro in `paper/numbers*.tex`, against the abstract text, and against the stored result files (including the d=9 certificate and its comparison with the recipes' intervals), with SHA-256 of the key stored files. Exit status is non-zero on any mismatch or missing file. `make verify-full` also re-runs `paper/make_numbers*.py` on the stored results in a temporary copy and requires all `numbers*.tex`, `tab_*.tex` and number ledgers to equal the arXiv v7 files.

This recomputes tables and counts from stored results; it does **not** repeat the simulations (stim/pymatching sampling, exact enumeration), which take CPU-days. The scripts that produced the stored results are in `code/` and `exact/`.

## Layout
- `code/` stim circuits and noise model (`truth.py`: rotated surface-code memory circuits, detector error model, fault-stratified importance-sampling reference estimator with the stopping rules), data generators (`fitdata.py`, `fitdata_bw.py`, drivers `run_cov.py`, `run_bw.py`), fitting recipes R1-R5 and coverage (`coverage.py`, `coverage_bw.py`, `analyze_main.py`, `analyze_bw.py`), the d=9 certified interval and betting bounds (`certified.py`, `cert_d9.py`, `cert_d9_vs_recipes.py`, `betting_monotone_check.py`), exact enumeration helpers, direct Monte Carlo (`direct.py`), controls A/B (surrogate and data-only), amendments (`amend_gemini.py`, `inflate_rolling.py`, `r4p_bw.py`, `zerofault_r2e.py`, `final_checks_v7.py`, `blind_forecasts.py`), structural A-0 checks (`a0_struct.py`, `a0_cost_probe.py`).
- `exact/` exact low-weight enumeration of the failure polynomial (`run_shells.py`, `assemble*.py`, manifests for d=3, 5, 7).
- `results/` all stored results: replicates, reference runs (`truth*`, `holdout_truth_d9*`), IS validation, exact-model certified validation (`cert_validation/`), d=9 certificate (`cert_d9_*`), coverage tables, controls, amendments. `logs/` run manifests and cost probes. `figs_amend/` amendment figures.
- `paper/` arXiv v7 source (`main.tex`, `supplement.tex`, tables, `numbers*.tex`, `refs.bib`, `figs/`), the two PDFs, `ABSTRACT.txt`, the builders `make_numbers*.py`, `make_tab_val.py`, figure scripts, and the number ledgers `numbers*.json` (each number with its source file and rule) and `CLAIMS.md` (generated number-to-source list).
- `plans/` the pre-registration texts (see below). `RESULTS*.md`, `REDTEAM_CHECKS_*.md` result notes.
- `CLAIMS.md`, `verify.py`, `Makefile`, `PATCHES.md`, `REDACTIONS.md`.

## Pre-registration plans and redactions
`plans/` holds `PLAN.md` as committed at each of its 26 commits (every amendment was committed before the runs it governs, as its heading states) and the frozen claims ledger, with the SHA-256 of each file in `plans/README.md`. In these texts and in the narrative notes, each internal operations note (machine names, launch, remote-access and monitoring instructions, screening remarks, local paths) is replaced by "[internal operations note redacted]". **The redactions are operational only; the originals are available to referees on request.** `REDACTIONS.md` lists the affected files. Scripts had absolute local paths and host names replaced (`PATCHES.md` gives the original and new SHA-256 of every patched file, because some files were pinned by hash prefix in the plans); no estimator, seed or statistic was edited.

## Licences
Code (`code/`, `exact/`, `paper/*.py`, `verify.py`): MIT, `LICENSE`. Paper text, figures, plans and stored results: CC BY 4.0, `LICENSE-DATA-TEXT`. No third-party data are used (all data are generated by the code in this bundle). Environment: Python 3.12, NumPy 2.5.5, SciPy 1.18.1 (`requirements.txt`); the simulations used stim 1.16.0 and pymatching 2.4.0.

## Known gaps
- `paper/main.tex` (v7) says in "Data and code availability" that the archive is not yet public and carries no URL; it needs the final repository URL/DOI before posting.
- `paper/tab_truth.tex` in the arXiv source has a hand-added z-score column that no script writes (see `PATCHES.md`); `--full` compares it without that column.
- Shell scripts that only launched or monitored jobs on specific machines were not included.
