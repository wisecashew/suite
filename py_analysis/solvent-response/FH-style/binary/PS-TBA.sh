#!/bin/bash

set -e

python nuanced_arr.py --label NECK --skips 1000 --draw-bin \
--xlim 0 0.4 --ylim 200 500 \
--Tbot 200 200 200 --Ttop 500 500 500 --T 200 500 --Tcrit 200 500 \
--vs 2.15333589 2.16532315 2.17391301 --vm 76.7661779 250.815736 599.786669 \
--pv 1 1 1 --pwms 0.63156 0.62208 0.614514 --pwmm 0 0 0 --pwss 0 0 0 \
--emsa -2496.22713 -2496.36601 -2496.31595 --emsn -100.119681 -100.090346 -100.077299 \
--emma -5006.82781 -5006.80995 -5006.84473 --emmn -5006.82781 -5006.80995 -5006.84473 \
--essa 0 0 0 --essn 0 0 0 \
--csv PS-tertbutylacetate_100k.csv PS-tertbutylacetate_233k.csv PS-tertbutylacetate_600k.csv \
--img PS-TBA.svg --dump PS-TBA-100k.dump PS-TBA-233k.dump PS-TBA-600k.dump 

:<<'END'
# this is the 100k good!
python constricted-cmaes-with_labels.py --label NECK --skip 1000 \
--T 10 500 --T-neck 350 --vm 76.7661779 --vs 0.919397956 \
--xfrmpv 23.0114718 --xfrmpw 0.538913872 \
--emsa -2496.22713 --emsn -100.119681 \
--emm -5006.82781 --ess 0 \
--img PS-TBA-cons-100k --csv PS-tertbutylacetate_100k.csv --maxiter 10 

# this is the 223k good!
python constricted-cmaes-with_labels.py --label NECK --skip 1000 \
--T 10 500 --T-neck 350 --vm 250.815736 --vs 0.889622589 \
--xfrmpv 23.0871052 --xfrmpw 0.498384478 \
--emsa -2496.36601 --emsn -100.090346 \
--emm -5006.80995 --ess 0 \
--img PS-TBA-cons-223k --csv PS-tertbutylacetate_233k.csv --maxiter 10 

# this is the 600k good!
python constricted-cmaes-with_labels.py --label NECK --skip 1000 \
--T 10 500 --T-neck 350 --vm 599.786669 --vs 0.846883104 \
--xfrmpv 23.0951364 --xfrmpw 0.466325210 \
--emsa -2496.31595 --emsn -100.077299 \
--emm -5006.84473 --ess 0 \
--img PS-TBA-cons-600k --csv PS-tertbutylacetate_600k.csv --maxiter 10 
END
