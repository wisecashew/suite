#!/bin/bash

set -e 

python binodal_parser.py --chisc -11 --chipc -1 --chips -1 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-1.0-chipc_-1.0-chisc_-11.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --island-stable-pkl IF_4/stable_islands/chips_-1-chipc_-1-chisc_-11-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-1-chipc_-1-chisc_-11-vs_1-vc_1-vp_32.unstable_islands.pkl --crit-pkl IF_4/crits/chips_-1-chipc_-1-chisc_-11-vs_1-vc_1-vp_32.crits.pkl --binodal-pkl IF_4/binodals/chips_-1.0-chipc_-1.0-chisc_-11.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --img chisc_dom

python binodal_parser.py --chisc -1 --chipc -11 --chips -1 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-1.0-chipc_-11.0-chisc_-1.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --island-stable-pkl IF_4/stable_islands/chips_-1-chipc_-11-chisc_-1-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-1-chipc_-11-chisc_-1-vs_1-vc_1-vp_32.unstable_islands.pkl --crit-pkl IF_4/crits/chips_-1-chipc_-11-chisc_-1-vs_1-vc_1-vp_32.crits.pkl --binodal-pkl IF_4/binodals/chips_-1.0-chipc_-11.0-chisc_-1.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --img chipc_dom

python binodal_parser.py --chisc -1 --chipc -1 --chips -11 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-11.0-chipc_-1.0-chisc_-1.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --island-stable-pkl IF_4/stable_islands/chips_-11-chipc_-1-chisc_-1-vs_1-vc_1-vp_32.stable_islands.pkl --island-unstable-pkl IF_4/unstable_islands/chips_-11-chipc_-1-chisc_-1-vs_1-vc_1-vp_32.unstable_islands.pkl --crit-pkl IF_4/crits/chips_-11-chipc_-1-chisc_-1-vs_1-vc_1-vp_32.crits.pkl --binodal-pkl IF_4/binodals/chips_-11.0-chipc_-1.0-chisc_-1.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --img chips_dom
