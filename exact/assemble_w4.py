"""Manifest for the SAMPLED weight-4 shell (see run_w4_sampled.py). usage: assemble_w4.py D
Reports (1) OBSERVED sums over the units done (hard lower bounds, counted failing sets only) and (2) a design-based ESTIMATE of the full-shell sums.
Units are processed in a seeded uniformly random order, so the done units are a simple random sample (without replacement) of the unit list; the
estimate is N_units * mean(unit sum) with the finite-population-corrected standard error. TRUTH GENERATION for the observed part; the estimate is an
estimate with a sampling error, NOT an exact coefficient."""
import sys, os, json, glob, hashlib, math, numpy as np
D = int(sys.argv[1]); OUT = os.path.expanduser(f"work/L_exact/out/d{D}")
des = json.load(open(f"{OUT}/w4s/design.json")); N = des["n_units"]
fs = sorted(glob.glob(f"{OUT}/w4s/*/u*.json")); us = [json.load(open(f)) for f in fs]
assert all(u["dem_sha"] == des["dem_sha"] and u["w"] == 4 for u in us), "unit/dem mismatch"
ids = [u["unit"] for u in us]; assert len(set(ids)) == len(ids)
n = len(us); sets = sum(u["n_sets"] for u in us); T = np.sum([u["T"] for u in us], axis=0); nf = sum(u["n_fail"] for u in us)
h = hashlib.sha256(b"".join(open(f, "rb").read() for f in fs)).hexdigest()
man = dict(kind="TRUTH GENERATION (exact decoding of every enumerated set); weight-4 shell PARTIAL and sampled, not complete", d=D, rounds=D, shell=4,
           design=des, units_done=n, units_total=N, frac_units_done=n / N, sets_done=sets, sets_total=des["total_sets"], frac_sets_done=sets / des["total_sets"],
           observed_lower_bounds=dict(n_fail=nf, T=T.tolist(), note="sums over processed units only; c_4 >= T[0] (given weights 1-3 have no failing sets)"), unit_files_sha256=h)
if n >= 30:
    est = {}
    for name, key in (("n_fail", lambda u: u["n_fail"]), ("T0", lambda u: u["T"][0]), ("T1", lambda u: u["T"][1]), ("T2", lambda u: u["T"][2]), ("T3", lambda u: u["T"][3])):
        x = np.array([key(u) for u in us], dtype=float); tot = N * x.mean(); se = N * x.std(ddof=1) / math.sqrt(n) * math.sqrt(max(1 - n / N, 0))
        est[name] = dict(estimate=tot, se=se, rel_se=se / tot if tot else None, ci95=[tot - 1.96 * se, tot + 1.96 * se])
    man["ESTIMATE_full_shell (design-based, simple random sample of units; NOT exact)"] = est
    man["coefficient_c4"] = dict(observed_lower_bound=T[0], estimate=est["T0"]["estimate"], se=est["T0"]["se"],
        status="ESTIMATE with sampling error (c_4 = sum over failing 4-sets of prod m_i, because weights 1-3 have zero failing sets)")
json.dump(man, open(f"{OUT}/manifest_w4.json", "w"), indent=1); print(json.dumps(man, indent=1))
