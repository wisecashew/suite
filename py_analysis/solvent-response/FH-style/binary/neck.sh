#!/bin/bash

set -e

# python neck_binodals.py --label NECK --skip 1000 \
# --xlim 0 0.4 --ylim 200 500 \
# --T 10 500 --T-neck 350 --vm 76.7661779 --vs 0.919397956 \
# --xfrmpv 23.0114718 --xfrmpw 0.538913872 \
# --emsa -2496.22713 --emsn -100.119681 \
# --emm -5006.82781 --ess 0 \
# --colors "coral" --marker "^" \
# --img PS-TBA-cons-100k --csv PS-tertbutylacetate_100k.csv 

python neck_binodals.py --label NECK --skip 100 \
--xlim 0 0.4 --ylim 200 500 \
--T 10 500 --T-neck 350 --vm 76.7661779 250.815736 599.786669 --vs 0.919397956 0.889622589 0.846883104 \
--xfrmpv 23.0114718 23.0871052 23.0951364 --xfrmpw 0.538913872 0.498384478 0.466325210 \
--emsa -2496.22713 -2496.36601 -2496.31595 --emsn -100.119681 -100.090346 -100.077299 \
--emm -5006.82781 -5006.80995 -5006.84473 --ess 0 0 0 \
--colors "coral" "red" "firebrick" --marker "^" "o" "s" \
--img PS-TBA.svg --csv PS-tertbutylacetate_100k.csv PS-tertbutylacetate_233k.csv PS-tertbutylacetate_600k.csv

:<<'END'
# this is the 223k good!
python neck_binodals.py --label NECK --skip 1000 \
--xlim 0 0.4 --ylim 200 500 \
--T 10 500 --T-neck 350 --vm 250.815736 --vs 0.889622589 \
--xfrmpv 23.0871052 --xfrmpw 0.498384478 \
--emsa -2496.36601 --emsn -100.090346 \
--emm -5006.80995 --ess 0 \
--img PS-TBA-cons-223k --csv PS-tertbutylacetate_233k.csv --maxiter 10 

# this is the 600k good!
python neck_binodals.py --label NECK --skip 1000 \
--xlim 0 0.4 --ylim 200 500 \
--T 10 500 --T-neck 350 --vm 599.786669 --vs 0.846883104 \
--xfrmpv 23.0951364 --xfrmpw 0.466325210 \
--emsa -2496.31595 --emsn -100.077299 \
--emm -5006.84473 --ess 0 \
--img PS-TBA-cons-600k --csv PS-tertbutylacetate_600k.csv --maxiter 10 
END
