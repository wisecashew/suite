#!/bin/bash

set -e
:<<'END'
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

python nuanced_arr.py --label LOOP --skips 500 --draw-bin \
--xlim 0 0.8 --ylim 300 600 \
--T 300 600 --Tbot 400 300 --Ttop 550 420 --Tcrit 400 550 \
--vs 5.37351454 1.74100733 --vm 167.315088 492.97102252 \
--pv 0.796872 1.0 --pwms 0.538365 0.265333 \
--pwmm 0 0 --pwss 0 0 \
--emsa -135.145959 -148.18038225 --emsn -54.7537934 -58.32471711 \
--emma -189.009097 -183.37638597 --emmn -189.009097 -183.37638597 \
--essa 0 0 --essn 0 0 \
--csv PEG-water.csv 2-butoxyethanol-water.csv --dump PEG-water.dump 2-butoxyethanol-water.dump \
--img loops.svg

# python constricted-cmaes-with_labels.py --label LOOP --T 10.0 1000.0 --T-loop 450.0 --skip 1000 --vm 167.03562589619366 --vs 7.918927221909584 --xfrmpv 0.6031562161928677 --xfrmpw 0.7491866417592598 --emm -186.52685910203877 --emsa -135.0444266897488 --emsn -53.8781130531603 --ess 0 --csv PEG-water.csv --img peg-water_over --sigma 0.01 --maxiter 200 # --optimize

# python constricted-cmaes-with_labels.py --label LOOP --T 10.0 800.0 --T-loop 350.0 --skip 1000 --vm 493.5459592772593 --vs 1.7664910252923334 --xfrmpv 27.842211646756688 --xfrmpw -1.0456883079773216 --emm -181.74463468310128 --emsa -147.61306886991608 --emsn -58.308171516843416 --ess 0 --csv 2-butoxyethanol-water.csv --img 2BE-water_over --sigma 0.01 --maxiter 200 # --optimize


