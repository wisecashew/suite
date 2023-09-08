#!/bin/bash

for chisc in -10 -7 -4 -1 0; do
    echo "chisc = $chisc..."
    python pchiexpl_v2.py -N 32 --chimeshsize 750 --nproc 96 --chiac ${chisc} --img cnscheck_chisc_${chisc}.png --llim -15 --ulim 15
done


