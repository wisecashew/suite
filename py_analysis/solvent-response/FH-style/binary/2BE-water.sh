#!/bin/bash

set -e

python nuanced_arr.py --label LOOP --skips 100 --draw-bin \
--xlim 0 0.3 --ylim 300 420 \
--Tbot 300 --Ttop 420 --T 300 420 --Tcrit 300 420 \
--vs 1.02231334 --vm 500.052965 \
--pv 1 \
--pwms 0.3521619336325539 \
--pwmm 0.3521619336325539 \
--pwss 0.3521619336325539 \
--emsa -533.725192 --emsn -101.667485 \
--emma -923.250513 --emmn -229.968308 \
--essa -0.0695408698 --essn 0.09234036 \
--csv 2-butoxyethanol-water.csv \
--img 2BE-water_cmaes_tot.svg > 2BE-water_cmaes.out 2>&1

