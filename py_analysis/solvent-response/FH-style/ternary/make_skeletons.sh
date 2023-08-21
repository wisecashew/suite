#!/bin/bash

set -e

for chibc in {1..12}; do
    echo "chibc = ${chibc}..."
    python pskelbin_v1.py --chiac 0 --chiab -1 --chibc -${chibc} --mesh 200 -N 32 --skelfile bin-N_32-chiab_-1-chibc_-${chibc}-chiac_0.skelfile
done

sbatch binodal_maker.slurm
