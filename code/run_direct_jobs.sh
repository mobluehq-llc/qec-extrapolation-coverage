#!/bin/bash
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1
cat > ../jobs_direct.txt <<J
7 1.5e-3 2007
7 2e-3 2008
7 3e-3 2009
7 4e-3 2010
5 5e-4 2021
5 1e-3 2022
5 2e-3 2023
5 3e-3 2024
5 4e-3 2025
3 1e-4 2031
3 2e-4 2032
3 5e-4 2033
3 1e-3 2034
3 2e-3 2035
3 4e-3 2036
J
cat ../jobs_direct.txt | xargs -P 8 -L 1 bash -c 'python3 direct.py $0 $1 $2 ../results/direct_d$0_p$1.json 1000 8e8 > ../logs/direct_d$0_p$1.log 2>&1'
echo ALLDONE > ../logs/direct_alldone
