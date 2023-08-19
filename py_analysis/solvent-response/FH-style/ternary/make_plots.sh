#!/bin/bash

set -e

python spinwbin_v4.py --chiac -11 --chiab -1 --chibc -1 --nproc 96 -N 32 --dumpfile improved.skelfile --bin-boundary binodal.out --image improved-binodal.png 
