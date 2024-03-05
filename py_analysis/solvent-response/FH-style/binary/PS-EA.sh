#!/bin/bash

set -e

# python nuanced.py --label LCST --skips 1000 --draw-spin --draw-bin --vs 1 --vm 170 --emsa -1002 --emsn -100 --emma -1700 --emmn -1000 --essa -300 --essn -200 --pv 1 --pwms 0.35 --pwmm 0.25 --pwss 0.4 --Tbot 400 --Ttop 500 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate-100k.csv --img PS-EA_100k.png

# this is okay
# python nuanced.py --label LCST --skips 1000 --draw-spin --vs 1 --vm 400 --emsa -1034 --emsn -100 --emma -1700 --emmn -1000 --essa -300 --essn -200 --pv 1 --pwms 0.35 --pwmm 0.35 --pwss 0.35 --Tbot 400 --Ttop 500 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate-100k.csv --img PS-EA_100k.png --draw-bin

# python nuanced.py --label LCST --skips 10 --draw-spin --vs 1 --vm 100 --emsa -1032 --emsn -100 --emma -1700 --emmn -1000 --essa -300 --essn -200 --pv 1 --pwms 0.35 --pwmm 0.35 --pwss 0.35 --Tbot 400 --Ttop 440 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate-100k.csv --img PS-EA_100k.png --draw-bin

# python nuanced.py --label LCST --skips 1000 --draw-spin --vs 1 --vm 250 --emsa -999.4 --emsn -100 --emma -1701 --emmn -1000 --essa -300 --essn -200 --pv 1 --pwms 0.35 --pwmm 0.25 --pwss 0.4 --Tbot 400 --Ttop 440 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate_223k.csv --img PS-EA_223k.png --draw-bin

# python nuanced.py --label LCST --skips 1000 --draw-spin --vs 1 --vm 300 --emsa -1028 --emsn -100 --emma -1700 --emmn -1000 --essa -300 --essn -200 --pv 1 --pwms 0.35 --pwmm 0.35 --pwss 0.35 --Tbot 400 --Ttop 440 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate_223k.csv --img PS-EA_223k.png --draw-bin


# python nuanced.py --label LCST --skips 1000 --draw-spin --vs 1 --vm 650 --emsa -1030.25 --emsn -100 --emma -1800 --emmn -1000 --essa -200 --essn -100 --pv 1 --pwms 0.35 --pwmm 0.35 --pwss 0.35 --Tbot 400 --Ttop 440 --T 400 500 --Tcrit 400 500 --csv PS-ethylacetate_600k.csv --img PS-EA_600k.png --draw-spin --draw-bin

python nuanced_arr.py --label LOOP --skips 10000 --draw-bin \
--xlim 0 0.3 --ylim 400 460 \
--Tbot 400 --Ttop 460 --T 400 450 --Tcrit 400 440 \
--vs 1.10606396 --vm 170.027238 \
--pv 1 \
--pwms 0.3530402658211769 \
--pwmm 0.3530402658211769 \
--pwss 0.3530402658211769 \
--emsa -799.9735 --emsn -131.521966   \
--emma -1289.99258 --emmn -169.924631 \
--essa -144.983301 --essn -50.0038813 \
--csv PS-ethylacetate-100k.csv \
--img PS-EA-100k-cmaes_tot.svg > PS-EA-100k-cmaes.out 2>&1



