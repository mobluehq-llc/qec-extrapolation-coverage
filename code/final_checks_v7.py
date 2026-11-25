"""Final-check (Gemini + ChatGPT, v7) stored-data checks. Stdlib only; run from this folder:  python3 code/final_checks_v7.py
Reads only stored results; writes results/L_final_checks_v7.json. One rule for every robust/undecided label, on INTEGER hit counts:
  k = hits of 300 at the reference point estimate; kmin, kmax = smallest / largest hits as the shared reference moves over +/-1.96 se.
  robust-below: kmax < 270 (270/300 = 90%); robust-at-or-above: kmin >= 270; undecided otherwise."""
import json, csv, math, statistics as st
R = "results"; N = 300; B90 = 270
TG = ["0.001", "0.0001", "1e-05", "1e-06"]
AM = json.load(open(f"{R}/L_amend_gemini.json"))["reference_uncertainty"]; BW = json.load(open(f"{R}/L_r4p_bw.json")); assert BW["_sanity_mismatch"] == 0
def cp(k, n=N, a=0.05):
    def cdf(x, p): return sum(math.comb(n, i) * p**i * (1-p)**(n-i) for i in range(0, x+1))
    def bis(f, lo, hi):
        for _ in range(80):
            mid = (lo+hi)/2
            if f(mid): hi = mid
            else: lo = mid
        return (lo+hi)/2
    lo = 0.0 if k == 0 else bis(lambda p: 1-cdf(k-1, p) >= a/2, 0, 1)
    hi = 1.0 if k == n else bis(lambda p: cdf(k, p) <= a/2, 0, 1)
    return lo, hi
def label(kmin, kmax): return "below" if kmax < B90 else ("above" if kmin >= B90 else "undecided")
cells = {}
for rc in ("R1", "R2", "R3", "R4p", "R5"):
    for d in ((7,) if rc == "R5" else (3, 5, 7)):
        for t in TG:
            key = f"{rc}|{d}|{t}|primary"; a = AM[key]; b = BW[key]
            for w, k, lo, hi in (("w1", round(a["hit_at_point"]*N), round(a["hit_min"]*N), round(a["hit_max"]*N)), ("w2", b["hits"], round(b["hit_min"]*N), round(b["hit_max"]*N))):
                assert abs(hi/N - (a["hit_max"] if w == "w1" else b["hit_max"])) < 1e-9
                l, u = cp(k); cells[f"{w}|{rc}|{d}|{t}"] = dict(k=k, kmin=lo, kmax=hi, cp=[l, u], label=label(lo, hi), point_below90=k < B90)
# frozen CSV cross-check of the CP interval
for w, f in (("w1", "L_coverage.csv"), ("w2", "L_bw_coverage.csv")):
    for r in csv.DictReader(open(f"{R}/{f}")):
        rc = {"R4p": "R4p"}.get(r["recipe"], r["recipe"]); key = f"{w}|{rc}|{r['d']}|{float(r['target']):g}"
        key = [x for x in cells if x.split("|")[:3] == key.split("|")[:3] and float(x.split("|")[3]) == float(r["target"])]
        if key and r["recipe"] in ("R1", "R2", "R3", "R5"):
            c = cells[key[0]]; assert c["k"] == int(r["hits"]), (key, c["k"], r["hits"]); assert abs(c["cp"][0]-float(r["cp_lo"])) < 1e-6 and abs(c["cp"][1]-float(r["cp_hi"])) < 1e-6, key
summ = {}
for w in ("w1", "w2"):
    for rc in ("R1", "R2", "R3", "R4p", "R5"):
        cs = {k: v for k, v in cells.items() if k.startswith(f"{w}|{rc}|")}
        for sub, ds in (("all", (3, 5, 7)), ("lowd", (3, 5)), ("d7", (7,))):
            sel = [v for k, v in cs.items() if int(k.split("|")[2]) in ds]
            summ[f"{w}|{rc}|{sub}"] = dict(n=len(sel), point_below90=sum(v["point_below90"] for v in sel), **{l: sum(v["label"] == l for v in sel) for l in ("below", "undecided", "above")})
# Table III rounding: how many printed cells put the point estimate outside its own INTEGER-rounded range
bad = [k for k, v in cells.items() if not (round(100*v["kmin"]/N) <= round(100*v["k"]/N, 1) <= round(100*v["kmax"]/N))]
# --- S7 caption medians (first window, scaled covariance) recomputed
def med(f, rcs):
    out = {}
    for rc in rcs:
        xs = [float(r["cov"])*100 for r in csv.DictReader(open(f"{R}/{f}")) if r["recipe"] == rc]; out[rc] = round(st.median(xs), 1)
    return out
s7 = dict(first_window=med("L_coverage_scaled.csv", ["R1", "R2", "R3", "R4", "R5"]), second_window=med("L_bw_coverage_scaled.csv", ["R1", "R2", "R3", "R4", "R5"]), caption="R1 8.0, R2 3.2, R3 51.0, R4 0.0, R5 15.7 (first window)")
# --- S14 sums
c9 = json.load(open(f"{R}/cert_d9_result.json"))["A0"]; per = c9["per"]
Ls = sum(p["P"]*p["lcb"] for p in per.values()); Us = sum(p["P"]*p["ucb"] for p in per.values()) + c9["tail"]
s14 = dict(L_sum=Ls, U_sum=Us, L_stored=c9["L"], U_stored=c9["U"], rows={w: dict(L=p["P"]*p["lcb"], U=p["P"]*p["ucb"]) for w, p in per.items()}, tail=c9["tail"], U_over_L=Us/Ls)
assert abs(Ls-c9["L"]) < 1e-20 and abs(Us-c9["U"]) < 1e-20
# --- interval scale: the 219% figure. Median log half-width h = 1.96*sigma; reported "half-width" = exp(h)-1 (upper); lower = 1-exp(-h)
h = math.log(1+2.19); scale = dict(upper_rel=math.exp(h)-1, lower_rel=1-math.exp(-h), U_over_L=math.exp(2*h), symmetric_halfwidth_UminusL_over_2p=(math.exp(h)-math.exp(-h))/2)
json.dump(dict(cells=cells, summary=summ, table3_integer_rounding_violations=bad, s7=s7, s14=s14, interval_scale=scale), open(f"{R}/L_final_checks_v7.json", "w"), indent=1)
print("S7", s7); print("S14 L,U", Ls, Us, Us/Ls); print("scale", scale); print("Table III cells with point outside integer-rounded range:", len(bad), bad)
for k, v in summ.items():
    if k.split("|")[1] in ("R2", "R3", "R4p"): print(k, v)
