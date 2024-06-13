#!/bin/bash

# for pv in 0.0 0.25 0.5 0.75 1.0; do
#         echo "pv = $pv..."
#         for g in 0 0.25 0.5 0.75 1.0; do
#                 echo -e "\tg = $g..."
#                 python heatmaps-vectorized-v3-medium-multicolor.py -g $g --pv $pv --image exploration-g_$g-pv_$pv.png &
#         done
#     wait
# done
pv=1.0
g=0.25
python heatmaps-vectorized-v3-singular-multicolor-deltas.py -g $g --pv $pv --image single-g_$g-pv_$pv.png --Emmn-llim -1.5 --Emmn-ulim 0 --Dmm-llim -0.5 --Dmm-ulim 0.5 --figsize 13.33 7.5
