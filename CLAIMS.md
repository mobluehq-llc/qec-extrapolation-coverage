# Claims ledger (L bundle, arXiv v7)
Headline numbers of the abstract (`paper/ABSTRACT.txt`, from `paper/main.tex`), each tied to (a) its macro in `paper/numbers*.tex`, (b) the abstract text, and (c) the stored result files.
`make verify` first **re-fits the registered recipes R1-R5 on the stored replicates** (`results/fit_d357.json`, `results/fit_bw.json`; 300 replicates x d=3,5,7 x five window points) in pure Python, requires the hit counts of all four coverage tables (52 cells each) to equal the stored CSVs, and recomputes the headline counts (24 of 24 and 19 of 24 cells below 90%; R2/R3 upper-limit excesses; R1 coverage on the second window) and compares them with the paper's macros. It then checks the rows below. `make verify-full` additionally regenerates every number and table with `paper/make_numbers*.py` (needs numpy and scipy).
Evidence paths are relative to this file; checks use the verify.py grammar. A check pattern uses `.` where the evidence line has a table bar.

| # | claim | evidence | check | status |
|---|---|---|---|---|
| 1 | Paper macro reps = 300 (replicates) | paper/numbers.tex | `grep:LV@reps\\endcsname\{300\}` | verified |
| 2 | Abstract states '300 replicates' | paper/ABSTRACT.txt | `grep:300\ replicates` | verified |
| 3 | Paper macro v6.fx.w1.pt = 24 (fixed-exponent cells below 90%, window 1) | paper/numbers_v6.tex | `grep:LV@v6\.fx\.w1\.pt\\endcsname\{24\}` | verified |
| 4 | Abstract states 'below 90% in 24 of 24 cells on the first window' | paper/ABSTRACT.txt | `grep:below\ 90%\ in\ 24\ of\ 24\ cells\ on\ the\ first\ window` | verified |
| 5 | Paper macro v6.fx.w2.pt = 19 (window 2) | paper/numbers_v6.tex | `grep:LV@v6\.fx\.w2\.pt\\endcsname\{19\}` | verified |
| 6 | Abstract states '19 of 24 on the second' | paper/ABSTRACT.txt | `grep:19\ of\ 24\ on\ the\ second` | verified |
| 7 | Paper macro v6.fx.w1.robust = 19 (below 90% throughout the reference range, window 1) | paper/numbers_v6.tex | `grep:LV@v6\.fx\.w1\.robust\\endcsname\{19\}` | verified |
| 8 | Abstract states 'in 19 and 11 cells' | paper/ABSTRACT.txt | `grep:in\ 19\ and\ 11\ cells` | verified |
| 9 | Paper macro v6.fx.w2.robust = 11 (window 2) | paper/numbers_v6.tex | `grep:LV@v6\.fx\.w2\.robust\\endcsname\{11\}` | verified |
| 10 | Abstract states 'in 19 and 11 cells' | paper/ABSTRACT.txt | `grep:in\ 19\ and\ 11\ cells` | verified |
| 11 | Paper macro hw.R2.max = 5.0 (R2 upper-limit excess max %) | paper/numbers.tex | `grep:LV@hw\.R2\.max\\endcsname\{5\.0\}` | verified |
| 12 | Abstract states 'R2, a constant per distance: 5.0--5.1%' | paper/ABSTRACT.txt | `grep:R2,\ a\ constant\ per\ distance:\ 5\.0\-\-5\.1%` | verified |
| 13 | Paper macro bw.hw.R2.max = 5.1 (R2, window 2) | paper/numbers_v4.tex | `grep:LV@bw\.hw\.R2\.max\\endcsname\{5\.1\}` | verified |
| 14 | Abstract states '5.0--5.1%' | paper/ABSTRACT.txt | `grep:5\.0\-\-5\.1%` | verified |
| 15 | Paper macro hw.R3.min = 2.8 (R3 min) | paper/numbers.tex | `grep:LV@hw\.R3\.min\\endcsname\{2\.8\}` | verified |
| 16 | Abstract states 'R3, one shared curve: 2.8--4.6%' | paper/ABSTRACT.txt | `grep:R3,\ one\ shared\ curve:\ 2\.8\-\-4\.6%` | verified |
| 17 | Paper macro hw.R3.max = 4.6 (R3 max) | paper/numbers.tex | `grep:LV@hw\.R3\.max\\endcsname\{4\.6\}` | verified |
| 18 | Abstract states '2.8--4.6%' | paper/ABSTRACT.txt | `grep:2\.8\-\-4\.6%` | verified |
| 19 | Paper macro bw.R1.covmin = 91.0 (R1 coverage min, window 2) | paper/numbers_v4.tex | `grep:LV@bw\.R1\.covmin\\endcsname\{91\.0\}` | verified |
| 20 | Abstract states 'covered 91.0--94.0% on the second' | paper/ABSTRACT.txt | `grep:covered\ 91\.0\-\-94\.0%\ on\ the\ second` | verified |
| 21 | Paper macro bw.R1.covmax = 94.0 (R1 coverage max, window 2) | paper/numbers_v4.tex | `grep:LV@bw\.R1\.covmax\\endcsname\{94\.0\}` | verified |
| 22 | Abstract states '91.0--94.0%' | paper/ABSTRACT.txt | `grep:91\.0\-\-94\.0%` | verified |
| 23 | Paper macro am.R4p.lowd.covmin = 94.7 (R4' coverage min d=3,5) | paper/numbers_amend.tex | `grep:LV@am\.R4p\.lowd\.covmin\\endcsname\{94\.7\}` | verified |
| 24 | Abstract states 'covered 94.7--95.7% at d=3 and 5' | paper/ABSTRACT.txt | `grep:covered\ 94\.7\-\-95\.7%\ at\ d=3\ and\ 5` | verified |
| 25 | Paper macro am.R4p.lowd.covmax = 95.7 (R4' coverage max) | paper/numbers_amend.tex | `grep:LV@am\.R4p\.lowd\.covmax\\endcsname\{95\.7\}` | verified |
| 26 | Abstract states '94.7--95.7%' | paper/ABSTRACT.txt | `grep:94\.7\-\-95\.7%` | verified |
| 27 | Paper macro ratio.fixed.lowpct.max = 27 (median prediction within % of reference, window 1) | paper/numbers.tex | `grep:LV@ratio\.fixed\.lowpct\.max\\endcsname\{27\}` | verified |
| 28 | Abstract states 'within 27% of the reference on the first window' | paper/ABSTRACT.txt | `grep:within\ 27%\ of\ the\ reference\ on\ the\ first\ window` | verified |
| 29 | Paper macro bw.fixed.maxdev = 21 (same, window 2) | paper/numbers_v4.tex | `grep:LV@bw\.fixed\.maxdev\\endcsname\{21\}` | verified |
| 30 | Abstract states '21% on the second' | paper/ABSTRACT.txt | `grep:21%\ on\ the\ second` | verified |
| 31 | Paper macro ratio.R1.d7.t6 = 80 (free-exponent fold error) | paper/numbers.tex | `grep:LV@ratio\.R1\.d7\.t6\\endcsname\{80\}` | verified |
| 32 | Abstract states 'high by up to 80-fold' | paper/ABSTRACT.txt | `grep:high\ by\ up\ to\ 80\-fold` | verified |
| 33 | Paper macro truthse.d7.far.min = 7 (d=7 reference error min %) | paper/numbers.tex | `grep:LV@truthse\.d7\.far\.min\\endcsname\{7\}` | verified |
| 34 | Abstract states 'reference errors 7--8% at p<=10^-4' | paper/ABSTRACT.txt | `grep:reference\ errors\ 7\-\-8%\ at\ p<=10\^\-4` | verified |
| 35 | Paper macro truthse.d7.far.max = 8 (max) | paper/numbers.tex | `grep:LV@truthse\.d7\.far\.max\\endcsname\{8\}` | verified |
| 36 | Abstract states 'reference errors 7--8%' | paper/ABSTRACT.txt | `grep:reference\ errors\ 7\-\-8%` | verified |
| 37 | Coverage table, first window: 52 cells (R1-R5, d=3,5,7, four targets) | results/L_coverage.csv | `count:^R[1-5], == 52` | verified |
| 38 | Coverage table, second window: 52 cells | results/L_bw_coverage.csv | `count:^R[1-5], == 52` | verified |
| 39 | First row of the first-window table (R1, d=3, p=1e-6): 42 of 300 replicate intervals contain the reference | results/L_coverage.csv | `csv:hits@recipe=R1 == 42` | verified |
| 40 | d=9 certificate: d = 9 | results/cert_d9_result.json | `json:d == 9` | verified |
| 41 | d=9 certificate: p = 1e-4 | results/cert_d9_result.json | `json:p == 0.0001` | verified |
| 42 | d=9 certificate (A-0 form): lower limit 6.52e-13 (3 s.f. rounded down: 6.52e-13) | results/cert_d9_result.json | `json:A0.L == 6.52e-13 +- 1e-15` | verified |
| 43 | d=9 certificate (A-0): upper limit 1.61e-11 | results/cert_d9_result.json | `json:A0.U == 1.61e-11 +- 1e-13` | verified |
| 44 | d=9 certificate: plug-in point estimate 3.8e-12 | results/cert_d9_result.json | `json:A0.est == 3.8e-12 +- 1e-13` | verified |
| 45 | d=9 certificate: assumption-free lower limit 6.30e-13 | results/cert_d9_result.json | `json:assumption_free.L == 6.30e-13 +- 1e-15` | verified |
| 46 | d=9 certificate: assumption-free upper limit 1.10e-6 | results/cert_d9_result.json | `json:assumption_free.U == 1.10e-6 +- 1e-8` | verified |
| 47 | d=9 certificate: 16 one-sided events in the union bound (A-0) | results/cert_d9_result.json | `json:A0.n_events == 16` | verified |
| 48 | d=9 certificate: 24 events, assumption-free | results/cert_d9_result.json | `json:assumption_free.n_events == 24` | verified |
| 49 | d=9 certificate: 37,126 failing draws in stratum 5 | results/cert_d9_result.json | `json:nfail.5 == 37126` | verified |
| 50 | d=9 certificate: 1e7 failure-free draws in each of strata 1-4 | results/cert_d9_result.json | `json:nsamples.1 == 10000000` | verified |
| 51 | R1 replicates whose whole 95% interval lies above the d=9 certified upper limit: 300 of 300 | results/cert_d9_vs_recipes.json | `json:A0.recipes.R1.certified_missed_high == 300` | verified |
| 52 | R2 replicates whose whole 95% interval lies above the d=9 certified upper limit: 0 of 300 | results/cert_d9_vs_recipes.json | `json:A0.recipes.R2.certified_missed_high == 0` | verified |
| 53 | R3 replicates whose whole 95% interval lies above the d=9 certified upper limit: 0 of 300 | results/cert_d9_vs_recipes.json | `json:A0.recipes.R3.certified_missed_high == 0` | verified |
| 54 | R4 replicates whose whole 95% interval lies above the d=9 certified upper limit: 239 of 300 | results/cert_d9_vs_recipes.json | `json:A0.recipes.R4.certified_missed_high == 239` | verified |
| 55 | R5: 208 of 300 above the certified upper limit | results/cert_d9_vs_recipes.json | `json:A0.recipes.R5.certified_missed_high == 208` | verified |
| 56 | Paper macro c9.L = 6.52\times10^{-13 | paper/numbers_v5.tex | `grep:LV@c9\.L\\endcsname\{6\.52\\times10\^\{\-13\}` | verified |
| 57 | Paper macro c9.U = 1.61\times10^{-11 | paper/numbers_v5.tex | `grep:LV@c9\.U\\endcsname\{1\.61\\times10\^\{\-11\}` | verified |
| 58 | Paper macro c9.est = 3.8\times10^{-12 | paper/numbers_v5.tex | `grep:LV@c9\.est\\endcsname\{3\.8\\times10\^\{\-12\}` | verified |
| 59 | Paper macro c9.free.U = 1.10\times10^{-6 | paper/numbers_v5.tex | `grep:LV@c9\.free\.U\\endcsname\{1\.10\\times10\^\{\-6\}` | verified |
| 60 | Stored-result integrity: results/fit_d357.json | results/fit_d357.json | `sha256:89f4f9322922c3557fc2e3cf8e4e5d3b61b1492ec01340808951dfbdf4b8a9f6` | verified |
| 61 | Stored-result integrity: results/fit_bw.json | results/fit_bw.json | `sha256:4db811282c8030a80b10978a22a34efcebe02d9dc7315aea20f4625fd4c733f0` | verified |
| 62 | Stored-result integrity: results/cert_d9_result.json | results/cert_d9_result.json | `sha256:6aae12c2bc79d0884fdcffdb4c2f466e4b6e0d8dcf130ad0a6448c3c6be57cb2` | verified |
| 63 | Stored-result integrity: results/cert_d9_vs_recipes.json | results/cert_d9_vs_recipes.json | `sha256:16fcfbd8243f024730a81fc0d6cb511d98c6f4e572ab987bbcf6e2fef2787cce` | verified |
| 64 | Stored-result integrity: results/L_amend_gemini.json | results/L_amend_gemini.json | `sha256:4a199ab05b22588ef2dd7f4a082f890f094265e9f789a254e519061d2ff9f56f` | verified |
| 65 | Stored-result integrity: results/L_coverage.csv | results/L_coverage.csv | `sha256:662da1cb501ad0a01f54dcf2294577bc1bc966d942d6c19b70945e14891a4715` | verified |
| 66 | Stored-result integrity: results/L_bw_coverage.csv | results/L_bw_coverage.csv | `sha256:f53c6e8e5ec1ad2473dd2a7008ee83a6cf4497616aaa9806c4909dc9f9df9644` | verified |
| 67 | Stored-result integrity: results/truth_d7_p1e-3.json | results/truth_d7_p1e-3.json | `sha256:37d518a19568910b10e2f1c78012389a227b7350cb22161681e42d0024c667c4` | verified |
