"""Assemble polynomial coefficients + manifest from unit files. usage: assemble.py D"""
import sys, os, json, glob, hashlib, math, numpy as np
from common import build, P_DEC
D = int(sys.argv[1]); OUT = os.path.expanduser(f"work/L_exact/out/d{D}")
EXTRA = sys.argv[2:]   # extra unit dirs (e.g. out_second/d5 pulled from the other machine); each unit i must appear exactly once overall
M = build(D); n = M["n"]; m = M["m"]
e = [1.0]  # elementary symmetric polys of all m, up to 6
poly = np.array([1.0]);
for x in m: poly = np.convolve(poly, [1.0, x])[:7]
e = poly.tolist()
W = {}; man = dict(kind="TRUTH GENERATION (exact enumeration), not a confirmatory test", d=D, rounds=D,
    n_mechanisms=n, n_detectors=M["nd"], model="independent DEM mechanisms q_i=m_i*p (m_i=q_i(1e-9)/1e-9, merged by stim, decompose_errors=True)",
    noise="after_clifford=before_measure_flip=after_reset_flip=p", decoder=f"pymatching, weights fixed from DEM at p={P_DEC}",
    dem_sha256=M["sha"], m_multiplicity_hist={str(k): int(v) for k, v in zip(*np.unique(np.round(m, 6), return_counts=True))}, shells={})
roots = [OUT] + [os.path.expanduser(x) for x in EXTRA]
ws_found = sorted({int(os.path.basename(x)[1:]) for r in roots for x in glob.glob(f"{r}/w*") if os.path.isdir(x) and os.path.basename(x)[1:].isdigit()})   # w4s (sampled shell) handled by assemble_w4.py
man["coverage"] = {}
for w in ws_found:
    seen = {}; dup = []
    for r in roots:
        for f in sorted(glob.glob(f"{r}/w{w}/u*.json")):
            i = int(os.path.basename(f)[1:6])
            if i in seen: dup.append(i)
            else: seen[i] = f
    missing = [i for i in range(n - w + 1) if i not in seen]
    man["coverage"][w] = dict(units_present=len(seen), units_total=n - w + 1, duplicates_ignored=dup, missing=len(missing), first_missing=missing[:5])
    assert not dup, f"w={w}: units on both machines {dup[:5]}"
    fs = [seen[i] for i in sorted(seen)]; us = [json.load(open(f)) for f in fs]
    assert all(u["dem_sha"] == M["sha"] and u["w"] == w and u["i"] == i for u, i in zip(us, sorted(seen))), "unit/dem mismatch"
    tot = math.comb(n, w); done = sum(u["n_sets"] for u in us)
    T = np.sum([u["T"] for u in us], axis=0); nf = sum(u["n_fail"] for u in us)
    h = hashlib.sha256(b"".join(open(f, "rb").read() for f in fs)).hexdigest()
    W[w] = (T, done == tot)
    man["shells"][w] = dict(units_done=len(fs), units_total=n - w + 1, sets_done=done, sets_total=tot, complete=done == tot,
        n_fail=nf, T=T.tolist(), unit_files_sha256=h, status="COMPLETE" if done == tot else f"PARTIAL ({done}/{tot} sets, lower bound, biased to low first-index units)")
coef = {}
for k in range(1, 7):
    ws = [w for w in W if w <= k]
    if not ws or any(w not in W for w in range(1, k + 1)): break
    c = 0.0
    for w in ws:
        r = k - w
        for j in range(r + 1): c += (-1) ** r * (-1) ** j * W[w][0][j] * e[r - j]
    ok = all(W[w][1] for w in ws)
    coef[k] = dict(coef=c, status="EXACT" if ok else "PARTIAL lower bound (counted failing sets only)")
man["coefficients_p_L=sum c_k p^k"] = coef
man["note"] = "c_k exact up to float64 summation; polynomial valid only where truncated tail is small (n_mech*p << 1)."
json.dump(man, open(f"{OUT}/manifest.json", "w"), indent=1); print(json.dumps(man, indent=1))
