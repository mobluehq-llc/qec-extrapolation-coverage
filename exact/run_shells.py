"""Exhaustive weight-shell enumeration with per-unit checkpoints (unit = first fault index).
usage: run_shells.py D WMAX [NPROC] [MAXHOURS]   Resumable: existing unit files are skipped."""
import sys, os, json, time, itertools, math, numpy as np
from multiprocessing import Pool
from common import build
D = int(sys.argv[1]); WMAX = int(sys.argv[2]); NP = int(sys.argv[3]) if len(sys.argv) > 3 else 12
DEADLINE = time.time() + 3600 * float(sys.argv[4]) if len(sys.argv) > 4 else 1e18
SHARD = os.environ.get("SHARD")   # "MOD:r1,r2,..." units with i%MOD in set belong to this machine (all shells >= MINW)
MINW = int(os.environ.get("MINW", "1")); SKIP = json.load(open(os.environ["SKIPFILE"])) if os.environ.get("SKIPFILE") else {}
if SHARD: _m, _r = SHARD.split(":"); SMOD = int(_m); SRES = {int(x) for x in _r.split(",")}
OUT = os.path.expanduser(f"work/L_exact/out/d{D}"); CH = 200_000
M = build(D); n = M["n"]; m = M["m"]; H = M["H"]; O = M["O"]; pm = M["pm"]

def unit_path(w, i): return f"{OUT}/w{w}/u{i:05d}.json"

def run_unit(a):
    w, i = a
    if time.time() > DEADLINE: return None
    nsets = math.comb(n - i - 1, w - 1); nfail = 0; T = np.zeros(4)
    k = w - 1
    it = itertools.combinations(range(i + 1, n), k) if k > 0 else iter([()])
    done = 0
    while done < nsets:
        take = min(CH, nsets - done)
        if k > 0:
            flat = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, take)), dtype=np.int32, count=take * k)
            idx = np.empty((take, w), dtype=np.int32); idx[:, 0] = i; idx[:, 1:] = flat.reshape(take, k)
        else:
            idx = np.full((1, 1), i, dtype=np.int32)
        S = H[idx[:, 0]].copy(); ob = O[idx[:, 0]].copy()
        for c in range(1, w): S ^= H[idx[:, c]]; ob ^= O[idx[:, c]]
        pred = pm.decode_batch(S).reshape(len(idx), -1)[:, 0]
        f = pred != ob
        if f.any():
            mm = m[idx[f]]; pr = mm.prod(1); p1 = mm.sum(1); p2 = (mm**2).sum(1); p3 = (mm**3).sum(1)
            h = [np.ones_like(p1), p1, (p1**2 + p2) / 2, (p1**3 + 3*p1*p2 + 2*p3) / 6]
            for j in range(4): T[j] += float((pr * h[j]).sum())
            nfail += int(f.sum())
        done += len(idx)
    r = dict(d=D, w=w, i=i, n_sets=nsets, n_fail=nfail, T=T.tolist(), dem_sha=M["sha"])
    os.makedirs(os.path.dirname(unit_path(w, i)), exist_ok=True)
    tmp = unit_path(w, i) + ".tmp"; json.dump(r, open(tmp, "w")); os.replace(tmp, unit_path(w, i))
    return (w, i, nsets, nfail)

if __name__ == "__main__":
    print(f"d={D} n_mech={n} n_det={M['nd']} sha={M['sha'][:12]}", flush=True)
    with Pool(NP) as pool:
        for w in range(1, WMAX + 1):
            todo = [(w, i) for i in range(0, n - w + 1) if not os.path.exists(unit_path(w, i))
                    and (not SHARD or w < MINW or i % SMOD in SRES) and i not in SKIP.get(str(w), [])]
            if w < MINW: todo = []
            t0 = time.time(); ns = 0
            print(f"shell w={w}: {len(todo)} units todo of {n-w+1}", flush=True)
            for r in pool.imap(run_unit, todo):   # ordered, largest unit first
                if r is None: print("deadline reached, stopping (resumable)", flush=True); sys.exit(0)
                ns += r[2]
                if r[1] % 10 == 0 or r[2] > 1e6:
                    print(f"  w={w} unit {r[1]} sets={r[2]} fail={r[3]} cum_sets={ns} {ns/(time.time()-t0):.0f}/s", flush=True)
            print(f"shell w={w} DONE {time.time()-t0:.0f}s", flush=True)
