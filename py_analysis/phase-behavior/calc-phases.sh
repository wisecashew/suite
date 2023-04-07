#!/bin/bash

for pv in 0.0 0.05 0.1 0.25 0.5 1.0; do
        echo "pv = $pv..."
        for g in 0 0.1 0.25 0.3 0.5 1.0; do
                echo -e "\tg = $g..."
                python heatmaps-cg.py -g $g -p $pv
        done
done
