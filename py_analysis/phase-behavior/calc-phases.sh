#!/bin/bash

for pv in 0.0 0.25 0.5 0.75 1.0; do
        echo "pv = $pv..."
        for g in 0 0.25 0.5 0.75 1.0; do
                echo -e "\tg = $g..."
                python heatmaps-vectorized-v3-medium-multicolor.py -g $g --pv $pv --image exploration-g_$g-pv_$pv.png &
        done
	wait
done
