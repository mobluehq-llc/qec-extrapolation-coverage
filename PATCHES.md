# Patches applied to the code for this bundle

Only absolute local paths, host names in comments and run instructions were changed (work directories now sit under `work/`; two driver scripts lost the host name from their file name and are now `code/run_cov.py` and `code/run_bw.py`). No estimator, seed, recipe or statistic was edited. Three further edits in `paper/make_numbers_v4.py` are described below.
Some files were pinned by SHA-256 prefix in `plans/` before their runs; those pins name the *original* files, so the original and bundle hashes are listed.

- `paper/make_numbers_v4.py`: a stale guard (`assert not exists(cert_d9_result.json)`, written while the d=9 run was in progress) was replaced by a comment, and the figure-only matplotlib import was made optional; its numbers and tables are unaffected.
- `paper/tab_truth.tex` in the arXiv v7 source carries a hand-added z-score column (keys `v8.z.*` from `make_numbers_v8.py`) that `make_numbers.py` does not write; `verify.py --full` compares it after removing that column.

| file | sha256 original | sha256 in this bundle |
|---|---|---|
| `paper/make_numbers_v3.py` | `8f787bda228dc845c003c51fe45893347321e9ad4e161cb6849843b8ef7c9bd2` | `cd5c80f2a0df5d044645762a92b5705ac4deea4b3a04493dbd03aff4fd8523b4` |
| `paper/make_figs.py` | `b3de899623c3cf6993b4dc73f47b4ff956fb91c37b744ca1a66bab2ee9663cb8` | `3904180b00323e96189c6e09384bfdaeea08971fc711650e12cf832a2768e744` |
| `paper/make_numbers.py` | `9e88fe05ceac28a4e863643f27be1330de1f19904f1fd18f73c2fa701ab53099` | `b2317d69c426d4a5c9bd557fe294232f6501408d131e636bdcc3279924548534` |
| `paper/make_numbers_v5.py` | `42c1ea0de141f558d195b1d92227332bcff695acc105161fce8a6b332c6786c3` | `355a82edca5098a6b9363f9da9a3546214a3b82a7cb449e714557ce8d4f9452b` |
| `paper/make_numbers_v4.py` | `4e8741ded7811515b352ad7f7663637a055c7ee6db083e324311716e6ea01df3` | `08e42e48d55b0fe9f2f27cb7577d4ec5e9b1515975785f5927aacfa79f41a20f` |
| `paper/make_numbers_amend.py` | `cdfc52233acc0b7786483478312758ecd2a77815efc36fa53b4d04d508393763` | `b03423d805e6f1c8002071a854262ecd1c63b25c8ce9da1f8fe3b08270539881` |
| `paper/make_fig_money_v3.py` | `4162c0882dc6a2e4b799c4e63fe68e147ff2549632f6160b3c3a0876bdcbe86a` | `74202e4408a4b5a35eb46db60949654083ef3c329cb5c1c7307ccbbcf340b6f6` |
| `code/redteam_checks_gemini.py` | `672df4a8ae2ae4826dacf3f5cace593e197d23bca0304e8e5cdc8a833a81f704` | `7ff9343b6705dc94e54bc8c4badc224b1ffc452679bb1bd27daf897a83f38729` |
| `code/controls.py` | `d1589a1bb1c8ce76ee2bda7d812693b9a451a3f94719bd2d87f72de8c0984344` | `329b2da598544ed0a06df8b3ad034a3744ca32ab9aac315e08c7350b26d7fffe` |
| `code/is_exact.py` | `a91c4e28132bad7a83d7026a3adcdb6bbf41816a06920f2ece52e6b6940759bb` | `ec4a8aa4df1ece1641fcf6f294c906da7c786d12b0985656629c014c8a8d718c` |
| `code/certified.py` | `32d97a577f829f1ffb588b6f2f4a42f4cad4c7894dfed5fb32beaeea995eb0eb` | `6cfa10c9ea5949dc91361c449eb2dad5bce57b469ee3a1f5bd72ec3d1a49edba` |
| `exact/run_w4_sampled.py` | `a47766040aad529cdfa545085a6ae0a5294364d7b9d92e2634ec322f7c2444c0` | `23b864a03c8a00259e9630fb62b6c544e401788e1a3b7b933afbd0842afa334e` |
| `exact/assemble_w4.py` | `cd88640ff33957fc18aba063c6539fe7b8c1676b44067f68c7993a73828690d3` | `066c7c16d8dabb73a6976071baf2a7cc0a8f71d8804c6eccc9f97ffd71a983f9` |
| `exact/assemble.py` | `0c5e72295d7682d006e4534eafc6154091bef97b05e5c61a18bff0cddcf64920` | `6f5cb70740aa5218b973e078a3f2a425979bd300572748ae47f7dd729ff1319f` |
| `exact/mc_validate.py` | `8daefd8ce2de6cc2883af951b45d21e6d1ec0e3c4298b5d12c07c5d227caf82d` | `8001fa3e277f0477e2a46c64eb8378fc4a70a80009a8758525ece3f2209d12e1` |
| `exact/run_shells.py` | `4fd4e5ed9d69531c4011f11b281259e884e998c582d3106adc86017616559130` | `92503beec73cefc9254885d0758428aef335b92b62326397430224d7ae14854a` |
| `code/betting_monotone_check.py` | `07cbd034bb6e090ab2dbfc95c5f0269f3f903e7e5a38abf20f6fa10710da4f7c` | `2ea56d667512233471aff52b3679784248266bdc2bf00eccfed237b64e49bcfb` |
| `code/zerofault_r2e.py` | `a9dcee59b106320e890b159e68f9cfe82e23b527cc0e8dc0d05f9515bbcc38b4` | `ae903776e0c6eff512b79f510c0c92a4728aff15b0abe16bd2df212b57402a01` |
| `code/run_cov.py` | `dca715984dd56cf97a0c5813daa26fab3a7e2f96f549c6371794bfdafb0e1139` | `50a4e0028ea55cbff9f708335bab0733f7624273c17bf962aa7edcd3f6f3d2c4` |
