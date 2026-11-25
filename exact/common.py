"""Shared model builder for L_exact. TRUTH GENERATION, not a confirmatory test.
Model: independent DEM mechanisms (stim DEM of surface_code:rotated_memory_z, rounds=d,
decompose_errors=True so mechanisms are merged by stim), mechanism i has probability q_i = m_i*p
with m_i the exact small-p multiplicity (q_i(1e-9)/1e-9). Decoder: pymatching, weights fixed from the
DEM at P_DEC (documented: a polynomial needs ONE fixed decoder)."""
import stim, numpy as np, pymatching, hashlib
P_LIN = 1e-9
P_DEC = 1e-3

def circuit(d, p):
    return stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=d,
        after_clifford_depolarization=p, before_measure_flip_probability=p,
        after_reset_flip_probability=p)

def mechs(dem):
    out = []
    for ins in dem.flattened():
        if ins.type != "error":
            continue
        dets = set(); obs = 0
        for t in ins.targets_copy():
            if t.is_relative_detector_id():
                dets ^= {t.val}
            elif t.is_logical_observable_id():
                obs ^= 1 << t.val
        out.append((ins.args_copy()[0], tuple(sorted(dets)), obs))
    return out

def build(d):
    dem_lin = circuit(d, P_LIN).detector_error_model(decompose_errors=True)
    dem_dec = circuit(d, P_DEC).detector_error_model(decompose_errors=True)
    ml, md = mechs(dem_lin), mechs(dem_dec)
    assert [(a, b) for _, a, b in ml] == [(a, b) for _, a, b in md], "DEM structure depends on p"
    n = len(ml); nd = dem_lin.num_detectors
    m = np.array([q / P_LIN for q, _, _ in ml])
    H = np.zeros((n, nd), dtype=np.uint8); O = np.zeros(n, dtype=np.uint8)
    for i, (_, dets, obs) in enumerate(ml):
        for j in dets: H[i, j] = 1
        O[i] = obs & 1
        assert obs in (0, 1)
    pm = pymatching.Matching.from_detector_error_model(dem_dec)
    sha = hashlib.sha256((str(dem_dec) + "|" + repr(list(m))).encode()).hexdigest()
    return dict(d=d, n=n, nd=nd, m=m, H=H, O=O, pm=pm, dem_dec=dem_dec, sha=sha)

def lin_dem(M, p):
    s = []
    for i in range(M["n"]):
        dets = np.nonzero(M["H"][i])[0]
        t = " ".join(f"D{j}" for j in dets) + (" L0" if M["O"][i] else "")
        s.append(f"error({float(M["m"][i]*p)!r}) {t}")
    return stim.DetectorErrorModel("\n".join(s))
