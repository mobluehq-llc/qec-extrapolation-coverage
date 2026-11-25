import json,glob,sys,math,csv,os
R=sys.argv[1]  # results dir containing truth_d*_p*.json and direct_d*_p*.json
rows=[]
for f in sorted(glob.glob(f'{R}/direct_d*_p*.json')):
    D=json.load(open(f)); t=f.replace('direct_','truth_')
    if not os.path.exists(t): continue
    T=json.load(open(t)); k,n=D['k'],D['n']; md=k/n; sdd=math.sqrt(max(k,1))/n
    z=(T['truth']-md)/math.hypot(T['se'],sdd)
    rows.append(dict(d=D['d'],p=D['p'],truth=T['truth'],truth_se=T['se'],direct=md,direct_se=sdd,direct_errors=k,shots=n,z=z,
                     eligible=int(T['truth']>=1e-6 and k>=100),ok=int(abs(z)<=2)))
w=csv.DictWriter(open(f'{R}/L_truth_validation.csv','w'),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
el=[r for r in rows if r['eligible']]
for r in rows: print(f"d={r['d']} p={r['p']:g} truth={r['truth']:.4e}+-{r['truth_se']:.1e} direct={r['direct']:.4e}+-{r['direct_se']:.1e} ({r['direct_errors']}) z={r['z']:+.2f} {'OK' if r['ok'] else 'FAIL'}")
print(f"VALIDATION: {sum(r['ok'] for r in el)}/{len(el)} eligible cells |z|<=2 = {sum(r['ok'] for r in el)/max(len(el),1):.3f}; gate {'PASS' if sum(r['ok'] for r in el)>=0.9*len(el) else 'FAIL'}")
