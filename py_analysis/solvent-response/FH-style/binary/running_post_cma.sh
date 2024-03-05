#!/bin/bash

set -e
:<<'END'
python nuanced.py --label UCST --draw-spin --draw-bin --skips 1000 \
--Tbot 275 --Ttop 350 --T 275 350 --Tcrit 275 350 \
--pv 0.402094 --pwms 0.41010049317 --pwmm 0.41010049317 --pwss 0.41010049317 \
--emsa -0.42745284 --emsn -1.15496282 \
--emma -13.768031 --emmn -13.31637158 \
--essa -0.40768503 --essn -0.22232703 \
--vm 151.09086636 --vs 1.25438629 \
--csv PIB-diisobutyl_ketone_227k.csv --img PIB-DIBK-227k-post-cma.png 

python nuanced.py --label UCST --draw-spin --draw-bin --skips 1000 \
--Tbot 275 --Ttop 350 --T 275 350 --Tcrit 275 350 \
--pv 0.431279 --pwms 0.6553688040962823 --pwmm 0.6553688040962823 --pwss 0.6553688040962823 \
--emsa -0.256374220 --emsn -0.205819368 \
--emma -14.2357021 --emmn -13.3154217 \
--essa -0.24642711 --essn 0.155786444 \
--vm 1000.01773 --vs 1.08880596 \
--csv PIB-diisobutyl_ketone_285k.csv --img PIB-DIBK-285k-post-cma.png 

python nuanced.py --label UCST --draw-spin --draw-bin --skips 1000 \
--Tbot 275 --Ttop 350 --T 275 350 --Tcrit 275 350 \
--pv 0.49323725203 --pwms 0.498172777004 --pwmm 0.498172777004 --pwss 0.498172777004 \
--emsa -0.0182113111 --emsn -0.0287151085 \
--emma -14.1778536 --emmn -14.2030666 \
--essa -0.00560155865 --essn 0.0392319384 \
--vm 4324.94932 --vs 0.997596506 \
--csv PIB-diisobutyl_ketone_6m.csv --img PIB-DIBK-6m-post-cma.png 

python nuanced.py --label LOOP --draw-spin --draw-bin --skips 1000 \
--Tbot 400 --Ttop 500 --T 400 500 --Tcrit 400 500 \
--pv 1 --pwms 0.352753 --pwmm 0.352753 --pwss 0.352753 \
--emsa -799.927351 --emsn -131.544062 \
--emma -1289.99867 --emmn -170.058725 \
--essa -144.994573 --essn -49.9946875 \
--vm 4324.94193 --vs 0.948831959 \
--csv PS-ethylacetate-100k.csv --img PS-ethylacetate-100k-post-cma.png 

python nuanced.py --label LOOP --draw-spin --draw-bin --skips 1000 \
--Tbot 400 --Ttop 550 --T 400 550 --Tcrit 400 550 \
--pv 0.999924795623999 --pwms 0.3682485455064879 --pwmm 0.3682485455064879 --pwss 0.3682485455064879 \
--emsa -720.195332 --emsn -141.935106 \
--emma -1167.32315 --emmn -130.927921 \
--essa -118.330748 --essn -13.8787480 \
--vm 156.088378 --vs 6.18112378 \
--csv PEG-water.csv --img PEG-water-post-cma.png
END

python nuanced.py --label NECK --draw-spin --draw-bin --skips 1000 \
--Tbot 200 --Ttop 500 --T 200 500 --Tcrit 200 500 \
--pv 1 --pwms 0.3534299393865181 --pwmm 0.3534299393865181 --pwss 0.3534299393865181 \
--emsa -2439.52880 --emsn -240.026102 \
--emma -4229.98574 --emmn -599.975512 \
--essa -768.011187 --essn -480.033377 \
--vm 75.0115804 --vs 0.998016369 \
--csv PS-tertbutylacetate_100k.csv --img PS-TBA-100k-post-cma.png

