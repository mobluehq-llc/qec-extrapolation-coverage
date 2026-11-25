"""MC of (a) linearised DEM model, (b) real stim circuit, decoded by the same fixed decoder.
usage: mc_validate.py D p shots [NPROC]  -> compares to polynomial from manifest."""
import sys, os, json, numpy as np, stim, math
from multiprocessing import Pool
from common import build, lin_dem, circuit
D = int(sys.argv[1]); P = float(sys.argv[2]); SH = int(float(sys.argv[3])); NP = int(sys.argv[4]) if len(sys.argv) > 4 else 12
M = build(D); pm = M["pm"]; B = 200_000
def work(a):
    kind, seed, shots = a
    s = (lin_dem(M, P).compile_sampler(seed=seed) if kind == "dem" else circuit(D, P).compile_detector_sampler(seed=seed))
    f = 0; done = 0
    while done < shots:
        k = min(B, shots - done)
        if kind == "dem": det, obs, _ = s.sample(k)
        else: det, obs = s.sample(k, separate_observables=True)
        pred = pm.decode_batch(det.astype(np.uint8)).reshape(k, -1)[:, 0]
        f += int((pred != obs[:, 0]).sum()); done += k
    return f
if __name__ == "__main__":
    man = json.load(open(os.path.expanduser(f"work/L_exact/out/d{D}/manifest.json")))
    c = man["coefficients_p_L=sum c_k p^k"]; poly = sum(v["coef"] * P**int(k) for k, v in c.items())
    res = dict(d=D, p=P, shots=SH, poly_value=poly, poly_orders={k: v["status"] for k, v in c.items()})
    for kind in ("dem", "circuit"):
        per = SH // (NP * 4); jobs = [(kind, 1000 * (1 + j) + (7 if kind == "dem" else 0), per) for j in range(NP * 4)]
        with Pool(NP) as pool: fl = sum(pool.map(work, jobs))
        N = per * len(jobs); ph = fl / N; se = math.sqrt(max(ph * (1 - ph), 1e-300) / N)
        res[kind] = dict(failures=fl, shots=N, p_hat=ph, se=se, z_vs_poly=(ph - poly) / se if se else None)
    print(json.dumps(res)); open(os.path.expanduser(f"work/L_exact/out/mc_d{D}_p{P}.json"), "w").write(json.dumps(res))
