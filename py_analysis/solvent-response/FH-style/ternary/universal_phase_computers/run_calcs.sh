#!/bin/bash

set -e

# python binodal_parser.py --final-binodal IF_2/chipc_emph/chips_-1.0-chipc_-9.0-chisc_-1.0-vs_5.0-vc_5.0-vp_2.0.binodals.pkl \
# --binodal-pkl IF_2/chipc_emph/chips_-1.0-chipc_-9.0-chisc_-1.0-vs_5.0-vc_5.0-vp_2.0.binodals.pkl \
# --chisc -1 --chipc -9 --chips -1 -vs 5 -vc 5 -vp 2 \
# --island-stable-pkl IF_2/stable_islands/chips_-1-chipc_-9-chisc_-1-vs_5-vc_5-vp_2.stable_islands.pkl \
# --island-unstable-pkl IF_2/unstable_islands/chips_-1-chipc_-9-chisc_-1-vs_5-vc_5-vp_2.unstable_islands.pkl \
# --crit-pkl IF_2/crits/chips_-1-chipc_-9-chisc_-1-vs_5-vc_5-vp_2.crits.pkl \
# --plot-edges --plot-crits --plot-binodals --img bin_test

# python binodal_parser.py --final-binodal IF_4/binodals/chips_-1.0-chipc_-10.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --binodal-pkl IF_4/binodals/chips_-1.0-chipc_-10.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --chips -1 --chipc -10 --chisc 0 -vs 1 -vc 1 -vp 32 --island-stable-pkl IF_4/stable_islands/chips_-1-chipc_-10-chisc_0-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-1-chipc_-10-chisc_0-vs_1-vc_1-vp_32.unstable_islands.pkl --plot-crits --crit-pkl IF_4/crits/chips_-1-chipc_-10-chisc_0-vs_1-vc_1-vp_32.crits.pkl --img chipc_-10-dom

# python binodal_parser.py --final-binodal IF_4/binodals/chips_-10.0-chipc_-1.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --binodal-pkl IF_4/binodals/chips_-10.0-chipc_-1.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --chips -10 --chipc -1 --chisc 0 -vs 1 -vc 1 -vp 32 --island-stable-pkl IF_4/stable_islands/chips_-10-chipc_-1-chisc_0-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-10-chipc_-1-chisc_0-vs_1-vc_1-vp_32.unstable_islands.pkl --plot-crits --crit-pkl IF_4/crits/chips_-10-chipc_-1-chisc_0-vs_1-vc_1-vp_32.crits.pkl --img chips_-10-dom

# python binodal_parser.py --final-binodal IF_4/binodals/chips_-1.0-chipc_-1.0-chisc_-10.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --binodal-pkl IF_4/binodals/chips_-1.0-chipc_-1.0-chisc_-10.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --chips -1 --chipc -1 --chisc -10 -vs 1 -vc 1 -vp 32 --island-stable-pkl IF_4/stable_islands/chips_-1-chipc_-1-chisc_10-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-1-chipc_-1-chisc_10-vs_1-vc_1-vp_32.unstable_islands.pkl --plot-crits --crit-pkl IF_4/crits/chips_-1-chipc_-1-chisc_-10-vs_1-vc_1-vp_32.crits.pkl --img chisc_-10-dom

chi=2.4
vs=1.5
vc=1.7
vp=1.2

python binodal_parser.py --final-binodal IF_1/binodals/chips_${chi}-chipc_${chi}-chisc_${chi}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl \
--binodal-pkl IF_1/binodals/chips_${chi}-chipc_${chi}-chisc_${chi}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl \
--island-stable-pkl IF_1/stable_islands/chips_${chi}-chipc_${chi}-chisc_${chi}-vs_${vs}-vc_${vc}-vp_${vp}.stable_islands.pkl \
--island-unstable-pkl IF_1/unstable_islands/chips_${chi}-chipc_${chi}-chisc_${chi}-vs_${vs}-vc_${vc}-vp_${vp}.unstable_islands.pkl \
--crit-pkl IF_1/crits/chips_${chi}-chipc_${chi}-chisc_${chi}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl --plot-edges \
--plot-crits --img diff_vols \
--chips ${chi} --chipc ${chi} --chisc ${chi} -vs ${vs} -vc ${vc} -vp ${vp}
