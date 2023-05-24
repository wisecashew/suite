#!/bin/bash

for chibc in -1 -2 -3 -4 -5 -6 -7 -8 -9 -10; do
        python ternary-plotter.py --chiab -1 --chibc ${chibc} --chiac 0 -N 32
done
