#!/bin/bash

for pv in 0.0 0.25 0.5 0.75 1.0; do
	echo "pv = $pv..."
	for g in 0 0.25 0.5 0.75 1.0; do
			echo -e "\tg = $g..."
			python heatmaps-vectorized-v3-medium-multicolor-deltas.py -g $g --pv $pv --image exploration-g_$g-pv_$pv.png --Emmn-llim -2 --Emmn-ulim 2 --Dmm-llim -2 --Dmm-ulim 2 &
	done
	wait
done

# for pv in 1.0; do
#         echo "pv = $pv..."
#        for g in 0.25; do
#                echo -e "\tg = $g..."
#                python heatmaps-vectorized-v3-medium-multicolor-deltas.py -g $g --pv $pv --image deltas-g_$g-pv_$pv.png --Emmn-llim -2 --Emmn-ulim 2 --Dmm-llim -2 --Dmm-ulim 2
#        done
#    wait
# done
