"""Weight-4 shell for d=7 (rounds 7): RANDOM-ORDER unit sampling. TRUTH GENERATION (exact decoding of every enumerated set), but the shell is NOT
completed: C(6602,4)=7.9e13 sets is far beyond ~24 h, so units are processed in a seeded random order and the result is labelled
PARTIAL: the observed sums are lower bounds, and assemble_w4.py also gives a design-based (Horvitz-Thompson) ESTIMATE of the full-shell sums, labelled ESTIMATE.
unit = (i, ja, jb): fault sets {i<j<k<l} with j in [ja,jb), sized to ~TARGET sets. Resumable (skips existing unit files), per-unit atomic checkpoints.
usage: run_w4_sampled.py D NPROC DEADLINE_EPOCH [TARGET]      (same model/decoder/sums as run_shells.py: T = sum over failing sets of prod(m) * h_j)"""
import sys, os, json, time, itertools, math, numpy as np
from multiprocessing import Pool
from common import build
D = int(sys.argv[1]); NP = int(sys.argv[2]); DEADLINE = float(sys.argv[3]); TARGET = float(sys.argv[4]) if len(sys.argv) > 4 else 2e8
SEED = 7; CH = 200_000
OUT = os.path.expanduser(f"work/L_exact/out/d{D}/w4s")
M = build(D); n = M["n"]; m = M["m"]; H = M["H"]; O = M["O"]; pm = M["pm"]

def make_units():
    per_j = np.array([math.comb(n - j - 1, 2) for j in range(n)], dtype=np.float64); units = []
    for i in range(n - 3):
        j = i + 1
        while j <= n - 3:
            ja = j; acc = 0.0
            while j <= n - 3 and (acc == 0.0 or acc + per_j[j] <= TARGET): acc += per_j[j]; j += 1
            units.append((i, ja, j))
    return units
UNITS = make_units()

def path(u): return f"{OUT}/{u // 1000:04d}/u{u:07d}.json"

def run_unit(u):
    if time.time() > DEADLINE: return None
    i, ja, jb = UNITS[u]; nsets = 0; nfail = 0; T = np.zeros(4)
    for j in range(ja, jb):
        S0 = H[i] ^ H[j]; o0 = O[i] ^ O[j]
        it = itertools.combinations(range(j + 1, n), 2); tot = math.comb(n - j - 1, 2); done = 0
        while done < tot:
            take = min(CH, tot - done)
            flat = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, take)), dtype=np.int32, count=take * 2).reshape(take, 2)
            S = S0[None, :] ^ H[flat[:, 0]] ^ H[flat[:, 1]]; ob = o0 ^ O[flat[:, 0]] ^ O[flat[:, 1]]
            pred = pm.decode_batch(S).reshape(take, -1)[:, 0]; f = pred != ob
            if f.any():
                fl = flat[f]; mm = np.concatenate([np.full((len(fl), 1), m[i]), np.full((len(fl), 1), m[j]), m[fl]], axis=1)
                pr = mm.prod(1); p1 = mm.sum(1); p2 = (mm**2).sum(1); p3 = (mm**3).sum(1)
                h = [np.ones_like(p1), p1, (p1**2 + p2) / 2, (p1**3 + 3*p1*p2 + 2*p3) / 6]
                for q in range(4): T[q] += float((pr * h[q]).sum())
                nfail += int(f.sum())
            done += take
        nsets += tot
    r = dict(d=D, w=4, unit=u, i=i, ja=ja, jb=jb, n_sets=nsets, n_fail=nfail, T=T.tolist(), dem_sha=M["sha"])
    os.makedirs(os.path.dirname(path(u)), exist_ok=True); tmp = path(u) + ".tmp"; json.dump(r, open(tmp, "w")); os.replace(tmp, path(u))
    return (u, nsets, nfail)

if __name__ == "__main__":
    order = np.random.default_rng(SEED).permutation(len(UNITS)).tolist()
    todo = [u for u in order if not os.path.exists(path(u))]
    print(f"d={D} n_mech={n} sha={M['sha'][:12]} units={len(UNITS)} (target {TARGET:.0e} sets) todo={len(todo)} seed={SEED}", flush=True)
    json.dump(dict(d=D, n=n, target=TARGET, seed=SEED, n_units=len(UNITS), total_sets=math.comb(n, 4), sets_in_units=int(sum(math.comb(n - j - 1, 2) for (i, ja, jb) in UNITS for j in range(ja, jb))), dem_sha=M["sha"]),
              open(f"{OUT}/design.json", "w") if os.makedirs(OUT, exist_ok=True) is None else None)
    t0 = time.time(); ns = 0; nu = 0
    with Pool(NP) as pool:
        for r in pool.imap_unordered(run_unit, todo, chunksize=1):
            if r is None: print("deadline reached, stopping (resumable)", flush=True); pool.terminate(); break
            ns += r[1]; nu += 1
            if nu % 20 == 0: print(f"  units_done_this_run={nu} sets={ns:.3e} {ns/(time.time()-t0):.0f}/s elapsed={time.time()-t0:.0f}s", flush=True)
