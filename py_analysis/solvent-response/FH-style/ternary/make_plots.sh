#!/bin/bash

for chiac in 1 1.25 1.5 1.75 2 2.25 2.5 3; do
        echo "Running calculations for chiac = $chiac."
        python plspin.py --chiac $chiac --chiab 1 --chibc 1 -N 32
done
