#!/bin/bash

set -e

# python binodal_parser.py --chisc 0 --chips -10 --chipc -1 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-10.0-chipc_-1.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --crit-pkl IF_4/crits/chips_-10-chipc_-1-chisc_0-vs_1-vc_1-vp_32.crits.pkl --plot-crits --img dom_chips.png

# python binodal_parser.py --chisc 0 --chips -1 --chipc -10 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-1.0-chipc_-10.0-chisc_0.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --crit-pkl IF_4/crits/chips_-1-chipc_-10-chisc_0-vs_1-vc_1-vp_32.crits.pkl --plot-crits --img dom_chipc.png

# python binodal_parser.py --chisc -10 --chips -1 --chipc -1 -vs 1 -vc 1 -vp 32 --final-binodal IF_4/binodals/chips_-1.0-chipc_-1.0-chisc_-10.0-vs_1.0-vc_1.0-vp_32.0.binodals.pkl --crit-pkl IF_4/crits/chips_-1-chipc_-1-chisc_-10-vs_1-vc_1-vp_32.crits.pkl --plot-crits --img dom_chisc.png

python cspin.py --chisc 0 --chips 1 --chipc 3 -vs 5 -vc 4 -vp 1 --ternary --plot-crits --tang_norm
python cspin.py --chisc 1 --chips 3 --chipc 1 -vs 1 -vc 3 -vp 1 --ternary --plot-crits --tang_norm
python cspin.py --chisc 2 --chips 1 --chipc 1 -vs 5 -vc 5 -vp 1 --ternary --plot-crits --tang_norm
python cspin.py --chisc 2 --chips 1 --chipc 1 -vs 5 -vc 5 -vp 1 --ternary --plot-crits --tang_norm
python cspin.py --chisc 2 --chips 3 --chipc 2 -vs 1 -vc 1 -vp 1 --ternary --plot-crits --tang_norm
python cspin.py --chisc 4 --chips 1 --chipc 1 -vs 1 -vc 1 -vp 5 --ternary --plot-crits --tang_norm

