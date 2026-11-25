#!/usr/bin/env python3
"""L bundle: (1) re-fit R1-R5 on the stored replicates in pure Python and recompute the coverage tables and headline counts; (2) re-run each CLAIMS.md row's check from its evidence file. Exit 1 on any mismatch,
unreadable evidence, or unverifiable row. Silence is never a pass (never silent).

CLAIMS.md row:  | # | claim | evidence | check | status |
  evidence = path (relative to CLAIMS.md) of the file the check reads
  check    = one of
     json:<dotted.path> <op> <value> [+-tol]    e.g.  json:coverage.e4r >= 0.95
                                                       json:fpr.E1 == 0.043 +- 0.0005
     csv:<column>@<col>=<val> <op> <value>      first row where col==val, read column
     grep:<regex>                               regex must match at least once in evidence
     count:<regex> == <n>                       exact number of matches
     sha256:<hex>                               hash of evidence file
     manual                                     FAILS (exit 1) unless --allow-manual; list is printed
  ops: == != >= <= > <
Status column: 'frozen'/'verified'/anything is ignored except 'dropped' (row skipped, and printed).
Usage: verify.py CLAIMS.md [--full] [--allow-manual]   (--full also regenerates every number/table; needs numpy + scipy)
"""
import sys, re, json, csv, hashlib, operator, math, glob, statistics
from pathlib import Path
OPS = {"==": operator.eq, "!=": operator.ne, ">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}

def rows(md):
    for ln, line in enumerate(md.splitlines(), 1):
        if not line.strip().startswith("|"): continue
        cells = [c.strip().strip("`") for c in line.strip().strip("|").split("|")]
        if len(cells) < 5 or set("".join(cells)) <= set("-: ") or cells[0] == "#": continue
        yield ln, cells[0], cells[1], cells[2].strip("`"), cells[3], cells[4]

def cmp(got, expr):
    m = re.match(r"\s*(==|!=|>=|<=|>|<)\s*(\S+)(?:\s*(?:\+-|±)\s*(\S+))?\s*$", expr)
    if not m: raise ValueError("bad comparison: " + expr)
    op, want, tol = m.group(1), m.group(2), m.group(3)
    try:
        g, w = float(got), float(want)
        if tol is not None and op == "==": return abs(g - w) <= float(tol), f"got {g} want {w} +-{tol}"
        return OPS[op](g, w), f"got {g} want {op} {w}"
    except (TypeError, ValueError):
        return OPS[op](str(got), want.strip("'\"")), f"got {got!r} want {op} {want}"

def run(path, check, base):
    ev = base / path
    if not ev.is_file(): return False, f"MISSING evidence file {path}"
    kind, _, rest = check.partition(":")
    if kind == "sha256":
        h = hashlib.sha256(ev.read_bytes()).hexdigest()
        ok = h == rest.strip().lower(); return ok, ("sha256 " + h[:16] + "..." if ok else f"MISMATCH want {rest.strip()} got {h}")
    text = ev.read_text(errors="replace")
    if kind == "grep":
        ok = re.search(rest, text, re.M) is not None
        return ok, "pattern found" if ok else f"MISMATCH pattern {rest!r} not found"
    if kind == "count":
        m = re.match(r"(.*)\s*==\s*(\d+)\s*$", rest)
        if not m: return False, "bad count check"
        n = len(re.findall(m.group(1).strip(), text, re.M))
        return n == int(m.group(2)), f"MISMATCH want {m.group(2)} got {n}" if n != int(m.group(2)) else f"count {n}"
    if kind == "json":
        m = re.match(r"([^\s=!<>]+)\s*(.*)$", rest)
        cur = json.loads(text)
        try:
            for p in m.group(1).split("."): cur = cur[int(p)] if isinstance(cur, list) else cur[p]
        except (KeyError, IndexError, ValueError, TypeError):
            return False, f"MISMATCH path {m.group(1)} not found"
        ok, msg = cmp(cur, m.group(2)); return ok, ("" if ok else "MISMATCH ") + msg
    if kind == "csv":
        m = re.match(r"([^@\s]+)@([^=\s]+)=(\S+)\s*(.*)$", rest)
        if not m: return False, "bad csv check"
        col, kc, kv, expr = m.groups()
        for r in csv.DictReader(text.splitlines()):
            if r.get(kc) == kv:
                if col not in r: return False, f"MISMATCH no column {col}"
                ok, msg = cmp(r[col], expr); return ok, ("" if ok else "MISMATCH ") + msg
        return False, f"MISMATCH no row with {kc}={kv}"
    return False, f"unknown check type {kind!r}"

# ---- pure-stdlib re-implementation of code/coverage.py and code/coverage_bw.py (fit_all, hit rule) used by verify.py
import json, math, glob, csv, statistics
TARGETS = [1e-3, 1e-4, 1e-5, 1e-6]
WINDOWS = {"w1": [4e-3, 6e-3, 8e-3, 1e-2, 1.5e-2], "bw": [1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3]}
def _inv(A):
    n = len(A); M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c])); M[c], M[p] = M[p], M[c]; pv = M[c][c]; M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c:
                f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]
def wls(X, y, w):
    n, k = len(X), len(X[0]); A = [[sum(X[i][a] * X[i][b] * w[i] for i in range(n)) for b in range(k)] for a in range(k)]; C = _inv(A)
    v = [sum(X[i][a] * w[i] * y[i] for i in range(n)) for a in range(k)]; b = [sum(C[a][c] * v[c] for c in range(k)) for a in range(k)]
    chi2 = sum(w[i] * (y[i] - sum(X[i][a] * b[a] for a in range(k))) ** 2 for i in range(n)); return b, C, chi2, n - k
def quad(g, C): return sum(g[a] * C[a][b] * g[b] for a in range(len(g)) for b in range(len(g)))
def fit_all(data, tp, PS, scaled=False):
    X0 = [math.log(p) for p in PS]; xs = math.log(tp); out = {}; Y = {}; Wt = {}
    for d, (k, n) in data.items():
        P = [a / b for a, b in zip(k, n)]; Y[d] = [math.log(p) for p in P]; Wt[d] = [a / (1 - p) for a, p in zip(k, P)]
    sc = lambda C, chi2, dof: [[v * max(1.0, chi2 / dof) for v in r] for r in C] if scaled and dof > 0 else C
    for d in data:
        y, w = Y[d], Wt[d]
        b, C, c2, dof = wls([[1, x] for x in X0], y, w); C = sc(C, c2, dof); g = [1, xs]; out[("R1", d)] = (sum(a * c for a, c in zip(g, b)), quad(g, C))
        h = (d + 1) / 2; b, C, c2, dof = wls([[1]] * 5, [yy - h * x for yy, x in zip(y, X0)], w); C = sc(C, c2, dof); out[("R2", d)] = (b[0] + h * xs, C[0][0])
        b, C, c2, dof = wls([[1, p / 1e-2, (p / 1e-2) ** 2] for p in PS], [yy - (d / 2) * x for yy, x in zip(y, X0)], w); C = sc(C, c2, dof); g = [1, tp / 1e-2, (tp / 1e-2) ** 2]
        out[("R4", d)] = (sum(a * c for a, c in zip(g, b)) + (d / 2) * xs, quad(g, C))
    Xs, ys, ws = [], [], []
    for d in data:
        h = (d + 1) / 2; Xs += [[1, h]] * 5; ys += [yy - h * x for yy, x in zip(Y[d], X0)]; ws += Wt[d]
    b, C, c2, dof = wls(Xs, ys, ws); C = sc(C, c2, dof)
    for d in data:
        h = (d + 1) / 2; g = [1, h]; out[("R3", d)] = (sum(a * c for a, c in zip(g, b)) + h * xs, quad(g, C))
    if 3 in data and 5 in data:
        m3, v3 = out[("R1", 3)]; m5, v5 = out[("R1", 5)]
        for d in data:
            if d >= 7: s = (d - 5) // 2; out[("R5", d)] = (m5 + s * (m5 - m3), (1 + s) ** 2 * v5 + s ** 2 * v3)
    return out
def load_truth(R):
    T = {}
    for f in glob.glob(f"{R}/truth_d*_p*.json"): x = json.load(open(f)); T[(x["d"], x["p"])] = x
    for f in glob.glob(f"{R}/truthextra_d*_p*.json"):
        y = json.load(open(f)); x = T[(y["d"], y["p"])]; w1, w2 = 1 / x["se"] ** 2, 1 / y["se"] ** 2
        x["truth"] = (w1 * x["truth"] + w2 * y["truth"]) / (w1 + w2); x["se"] = (w1 + w2) ** -0.5
    return {k: v["truth"] for k, v in T.items()}
def load_reps(path):
    reps = {}
    for x in json.load(open(path)):
        r = reps.setdefault(x["rep"], {}).setdefault(x["d"], ([0] * 5, [0] * 5)); r[0][x["ip"]] = x["k"]; r[1][x["ip"]] = x["n"]
    return [reps[r] for r in sorted(reps)]
def coverage(R, window, scaled):
    PS = WINDOWS[window]; truth = load_truth(R); reps = load_reps(f"{R}/" + ("fit_d357.json" if window == "w1" else "fit_bw.json"))
    cnt, tot, hw = {}, {}, {}
    for data in reps:
        for tp in TARGETS:
            for (rc, d), (mu, var) in fit_all({d: data[d] for d in (3, 5, 7)}, tp, PS, scaled).items():
                t = truth.get((d, tp))
                if t is None: continue
                key = (rc, d, tp); cnt[key] = cnt.get(key, 0) + int(abs(math.log(t) - mu) <= 1.96 * math.sqrt(var)); tot[key] = tot.get(key, 0) + 1
                hw.setdefault(key, []).append(math.exp(1.96 * math.sqrt(var)) - 1)
    return cnt, tot, hw


def macro(base, key):
    for f in sorted((base / "paper").glob("numbers*.tex")):
        m = re.search(r"LV@" + re.escape(key) + r"\\endcsname\{([^}]*(?:\{[^}]*\}[^}]*)*)\}", f.read_text())
        if m: return m.group(1)
    raise KeyError(key)

def recompute(base):
    """(1) re-fit R1-R5 on the stored replicates in pure Python (no numpy) and require the per-cell hit counts of all four coverage tables to be identical to the stored CSVs;
    (2) recompute the headline counts and compare them with the paper's macros."""
    R = str(base / "results"); bad = 0; res = {}
    for w, s, fn in (("w1", False, "L_coverage.csv"), ("w1", True, "L_coverage_scaled.csv"), ("bw", False, "L_bw_coverage.csv"), ("bw", True, "L_bw_coverage_scaled.csv")):
        cnt, tot, hw = coverage(R, w, s)
        rows = {(r["recipe"], int(r["d"]), float(r["target"])): r for r in csv.DictReader(open(f"{R}/{fn}"))}
        mism = [k for k in rows if cnt.get(k) != int(rows[k]["hits"]) or tot.get(k) != int(rows[k]["reps"])]
        miscal = [k for k in rows if (cnt[k] / tot[k] < 0.90) != bool(int(rows[k]["miscal"]))]
        ok = not mism and not miscal and len(cnt) == len(rows) == 52
        print(("PASS  " if ok else "FAIL  ") + f"recompute: {fn}: {len(rows)} cells, hit counts from {tot[next(iter(tot))]} stored replicates identical" + ("" if ok else f" (MISMATCH {mism[:3]} {miscal[:3]})")); bad += (not ok)
        if not s: res[w] = (cnt, tot, hw)
    def chk(label, got, key, fmt):
        nonlocal bad
        want = macro(base, key); ok = fmt(got) == want; bad += (not ok)
        print(("PASS  " if ok else "FAIL  ") + f"recompute: {label} = {fmt(got)} (paper macro {key} = {want})" + ("" if ok else " MISMATCH"))
    for w, tag in (("w1", "first window"), ("bw", "second window")):
        cnt, tot, hw = res[w]; fx = [k for k in cnt if k[0] in ("R2", "R3")]
        chk(f"{tag}: fixed-exponent cells (R2, R3; 24) with coverage below 90% at the reference point", sum(cnt[k] / tot[k] < 0.9 for k in fx), "v6.fx.%s.pt" % w.replace("bw", "w2").replace("w1", "w1"), str)
        med = lambda k: statistics.median(hw[k]) * 100
        chk(f"{tag}: largest median upper-limit excess of R2", max(med(k) for k in hw if k[0] == "R2"), "hw.R2.max" if w == "w1" else "bw.hw.R2.max", lambda v: "%.1f" % v)
    cnt, tot, hw = res["w1"]; chk("first window: smallest median upper-limit excess of R3", min(statistics.median(hw[k]) * 100 for k in hw if k[0] == "R3"), "hw.R3.min", lambda v: "%.1f" % v)
    chk("first window: largest median upper-limit excess of R3", max(statistics.median(hw[k]) * 100 for k in hw if k[0] == "R3"), "hw.R3.max", lambda v: "%.1f" % v)
    cnt, tot, hw = res["bw"]; r1 = [cnt[k] / tot[k] * 100 for k in cnt if k[0] == "R1"]
    chk("second window: lowest R1 coverage (%)", min(r1), "bw.R1.covmin", lambda v: "%.1f" % v); chk("second window: highest R1 coverage (%)", max(r1), "bw.R1.covmax", lambda v: "%.1f" % v)
    return bad

def regen(base):
    """--full: re-run every paper/make_numbers*.py on the stored results in a temp copy of the bundle (needs numpy and scipy) and require numbers*.tex and tab_*.tex to equal the shipped arXiv v7 files."""
    import subprocess, tempfile, shutil, os
    try: import numpy, scipy
    except ImportError: print("FAIL  --full needs numpy and scipy (pip install -r requirements.txt)"); return 1
    tmp = Path(tempfile.mkdtemp()) / "b"; shutil.copytree(base, tmp, ignore=shutil.ignore_patterns("work")); bad = 0
    for s in ("make_numbers", "make_numbers_amend", "make_numbers_v3", "make_numbers_v4", "make_numbers_v5", "make_numbers_v6", "make_numbers_v7", "make_numbers_v8", "make_tab_val"):
        r = subprocess.run([sys.executable, "-I", s + ".py"], cwd=tmp / "paper", capture_output=True, text=True)
        if r.returncode: print("FAIL  regen:", s, "exited", r.returncode, r.stderr[-300:]); bad += 1
    norm = lambda p: re.sub(r"; \\V\{v8\.z\.[^}]*\}", "", p.read_text())
    names = sorted(p.name for p in (base / "paper").glob("numbers*.tex")) + sorted(p.name for p in (base / "paper").glob("tab_*.tex"))
    for nme in names:
        if norm(tmp / "paper" / nme) != norm(base / "paper" / nme): print("FAIL  regen: " + nme + " differs from the shipped arXiv v7 file"); bad += 1
    for f in sorted((base / "paper").glob("numbers*.json")):
        if json.loads((tmp / "paper" / f.name).read_text()) != json.loads(f.read_text()): print("FAIL  regen: " + f.name + " differs"); bad += 1
    if not bad: print(f"PASS  regen: {len(names)} numbers/table files and {len(list((base/'paper').glob('numbers*.json')))} number ledgers regenerated from stored results, identical to the arXiv v7 files")
    shutil.rmtree(tmp.parent); return bad

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    allow_manual = "--allow-manual" in sys.argv
    if not args: sys.exit("usage: verify.py CLAIMS.md [--allow-manual]")
    cp = Path(args[0]); base = cp.parent
    n = 1; bad = recompute(base); manual = []
    if "--full" in sys.argv: n += 1; bad += (regen(base) > 0)
    for ln, num, claim, ev, check, status in rows(cp.read_text()):
        if status.lower() == "dropped": print(f"SKIP  #{num} (dropped)"); continue
        n += 1
        if check.strip() == "manual" or not check.strip():
            manual.append(num); print(f"MANUAL #{num}: {claim[:60]}")
            if not allow_manual: bad += 1
            continue
        try: ok, msg = run(ev, check.strip(), base)
        except Exception as e: ok, msg = False, f"ERROR {type(e).__name__}: {e}"
        print(("PASS  " if ok else "FAIL  ") + f"#{num} [{ev}] {msg}")
        bad += (not ok)
    if n == 0: print("FAIL: no checkable rows in CLAIMS.md (cannot measure)"); sys.exit(1)
    print(f"\n{n - bad}/{n} rows verified; {bad} failing; manual rows: {manual or 'none'}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__": main()
