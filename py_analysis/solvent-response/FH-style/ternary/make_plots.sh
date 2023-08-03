#!/bin/bash

set -e

for chiac in 1 1.25 1.5 1.75 2 2.25 3 5; do
    echo "Running calculations for chibc = $chibc, chiab = $chiab."
    python plspin.py --chiac $chiac --chiab 1 --chibc 1 -N 32
done
