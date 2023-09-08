#!/bin/bash

set -e

for chisc in 11; do
    python pskelbin_uni1.py --chisc -${chisc} --chips -1 --chipc -1 --mesh 300 -N 32 --skelfile bin-N_32-chips_-1-chipc_-1-chisc_-${chisc}.skelfile > skel_for_-${chisc}.out 2>&1 &
done 
wait

for chisc in 11; do
    python spinwbin_g2critpoint.py --chisc -${chisc} --chips -1 --chipc -1 -N 32 --dumpfile bin-N_32-chips_-1-chipc_-1-chisc_-${chisc}.skelfile --bin-boundary N_32-chips_-1-chipc_-1-chisc_-${chisc}.binodal --tielines --image plot-N_32-chips_-1-chipc_-1-chisc_-${chisc}.png --ternary --tieline-density 100 > bin_for_-${chisc}.out 2>&1 & 
done
wait
