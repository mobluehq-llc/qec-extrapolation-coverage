Exact-truth generation (TRUTH GENERATION, not a confirmatory test). [internal operations note redacted]
Run: `python3 run_shells.py D WMAX 12 HOURS` (resumable, per-unit checkpoints in out/dD/wW/), `python3 assemble.py D` (manifest.json with counts + sha256), `python3 mc_validate.py D p shots`.
Model: independent stim DEM mechanisms q_i=m_i*p (m_i = small-p multiplicity), fixed pymatching decoder weighted at p=1e-3.

Two-machine split (d=5 w4): SHARD=22:r,... MINW=4 SKIPFILE=skip_w4.json; [internal operations note redacted] Merge: [internal operations note redacted]

## d=7 (rounds 7) run, started 2026-11-23 (TRUTH GENERATION; weights 1-3 target complete, weight 4 PARTIAL/sampled)
n_mechanisms=6602, n_detectors=336, dem_sha256 prefix be3a85cce1b8. C(n,w): w=1 6.6e3, w=2 2.2e7, w=3 4.8e10, w=4 7.9e13 (w=4 cannot be completed: ~1e3 times the 24 h budget).
- Weights 1-3: `python3 run_shells.py 7 3 14 24` (14 procs, resumable; unit files out/d7/w{1,2,3}/). Weights 1 and 2 finished with 0 failing sets.
- Weight 4: `python3 run_w4_sampled.py 7 14 <deadline_epoch> 2e8` (chained after w1-3 by `chain_w4.sh`, hard stop 2026-11-24 22:50). Units (i, j-block) of ~2e8 sets (413,8xx units) processed in a seeded RANDOM order, so the finished units are a simple random sample; unit files out/d7/w4s/. NOT the same unit scheme as run_shells.py, hence a separate script.
  Validated on d=3 w=4 (full shell): the sampled runner, run to completion with resume in the middle, gives sets 272,963,405, failing 95,856,732, T0 3623956.158895363 = the run_shells.py manifest (3623956.158895362).
- Assemble: `python3 assemble.py 7` (shells 1-3 -> out/d7/manifest.json) and `python3 assemble_w4.py 7` (out/d7/manifest_w4.json: OBSERVED lower bounds with unit-file sha256, plus a labelled design-based ESTIMATE of the full-shell sums with finite-population SE). The estimate is not an exact coefficient; c_4 is exact only if the shell were complete, which it is not.
- Progress: [internal operations note redacted]
