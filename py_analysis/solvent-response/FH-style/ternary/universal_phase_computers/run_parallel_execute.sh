#!/bin/bash

set -e 

python parallel_execute.py \
--addr-stable IF_7/REGULAR/stable_islands/. \
--addr-unstable IF_7/REGULAR/unstable_islands/. \
--addr-binodals IF_7/REGULAR/binodals/. --nproc 96 \
--addr-crits IF_7/REGULAR/crits/. --mesh mesh.pkl \
--binodal-pkl IF_7/REGULAR/binodals/.
