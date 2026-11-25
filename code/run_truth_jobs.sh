#!/bin/bash
# truth cells, 6 parallel. job = "d p seed"
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat > ../jobs_truth.txt <<J
7 1e-3 1007
7 1e-4 1008
7 1e-5 1009
7 1e-6 1010
7 4e-3 1011
7 3e-3 1012
7 2e-3 1013
7 1.5e-3 1014
5 1e-3 1021
5 1e-4 1022
5 1e-5 1023
5 1e-6 1024
5 5e-4 1025
5 2e-3 1026
5 3e-3 1027
5 4e-3 1028
3 1e-3 1031
3 1e-4 1032
3 1e-5 1033
3 1e-6 1034
3 2e-4 1035
3 5e-4 1036
3 2e-3 1037
3 4e-3 1038
J
cat ../jobs_truth.txt | xargs -P 6 -L 1 bash -c 'REL=0.02 python3 truth.py $0 $1 $2 ../results/truth_d$0_p$1.json 2400 > ../logs/truth_d$0_p$1.log 2>&1'
echo ALLDONE > ../logs/truth_alldone
