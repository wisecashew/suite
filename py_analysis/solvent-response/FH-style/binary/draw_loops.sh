#!/bin/bash

set -e
python nuanced_arr.py --label LOOP --skips 10000 --draw-bin \
--xlim 0 0.6 --ylim 300 550 \
--Tbot 400 --Ttop 550 --T 400 550 --Tcrit 400 550 \
--vs 6.18818954 --vm 156.149842 \
--pv 0.999925168910469 \
--pwms 0.3682443914686681 \
--pwmm 0.3682443914686681 \
--pwss 0.3682443914686681 \
--emsa -720.194421 --emsn -141.971716 \
--emma -1167.2924 --emmn -130.962698 \
--essa -118.331189 --essn -13.9077579 \
--csv PEG-water.csv --dump PEG-water.dump \
--img peg-water.png > peg-water.out 2>&1

:<<'END'
python nuanced_arr.py --label LOOP --skips 500 --draw-bin \
--xlim 0 0.7 --ylim 300 550 \
--Tbot 400 300 --Ttop 550 420 --T 300 550 --Tcrit 400 550 \
--vs 6.09750925 1.02725865 --vm 156.124071 500.016557 \
--pv 0.999922 1 \
--pwms 0.367884 0.351223 \
--pwmm 0.367884 0.351223 \
--pwss 0.367884 0.351223 \
--emsa -720.210927 -533.73979 --emsn -141.873416 -101.796721 \
--emma -1167.24828 -922.889237 --emmn -130.93181 -229.992313 \
--essa -118.30791 -0.115343225 --essn -13.8431102 0.092407127 \
--csv PEG-water.csv 2-butoxyethanol-water.csv \
--img loops.svg
END

