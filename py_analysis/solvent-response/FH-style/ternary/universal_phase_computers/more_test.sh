#!/bin/bash

set -e

python bcalculator.py \
--chisc 4 --chips 0 --chipc 1 -vs 4 -vc 4 -vp 1 \
--crit-pkl IF_7/crits/chisc_4_chips_0_chipc_1-vs_4-vc_4-vp_1.crits.pkl --mesh mesh.pkl \
--final-binodal my_bin.pkl --crit-pkl my_crit.pkl --img my_img --plot-crits --plot-binodal \
--island-stable-pkl IF_7/stable_islands/chisc_4_chips_0_chipc_1-vs_4-vc_4-vp_1.stable_islands.pkl \
--island-unstable-pkl IF_7/unstable_islands/chisc_4_chips_0_chipc_1-vs_4-vc_4-vp_1.unstable_islands.pkl 
