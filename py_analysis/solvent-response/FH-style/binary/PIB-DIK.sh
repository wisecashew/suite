#!/bin/bash

set -e

# python nuanced.py --label UCST --skips 1000 --draw-spin --draw-bin --vs 1 --vm 4325 --emsa 0 --emsn 0 --emma -14.2 --emmn -14.2 --essa 0 --essn 0 --pv 0 --pwms 0 --pwmm 0 --pwss 0 --Tbot 300 --Ttop 400 --T 300 400 --Tcrit 300 400 --csv PIB-diisobutyl_ketone_6m.csv --img PIB-DIBK_6m.png
# python nuanced.py --label UCST --skips 1000 --draw-spin --draw-bin --vs 1 --vm 1000 --emsa 0 --emsn 0 --emma -14.2 --emmn -14.2 --essa 0 --essn 0 --pv 0 --pwms 0 --pwmm 0 --pwss 0 --Tbot 300 --Ttop 400 --T 300 400 --Tcrit 300 400 --csv PIB-diisobutyl_ketone_285k.csv --img PIB-DIBK_285k.png
# python nuanced.py --label UCST --skips 1000 --draw-spin --draw-bin --vs 1 --vm 150 --emsa 0 --emsn 0 --emma -14.25 --emmn -14.25 --essa 0 --essn 0 --pv 0 --pwms 0 --pwmm 0 --pwss 0 --Tbot 250 --Ttop 350 --T 250 350 --Tcrit 250 350 --csv PIB-diisobutyl_ketone_227k.csv --img PIB-DIBK_227k.png

python nuanced_arr.py --label UCST --skips 10000 --draw-bin \
--xlim 0 0.4 --ylim 280 340 \
--Tbot 200 200 200 --Ttop 400 400 400 --T 200 400 --Tcrit 275 350 \
--vs 1.25438629 1.09213698 0.997596506 --vm 151.09086636 1000.01225 4324.94706 \
--pv 0.4020938049827699 0.4312796638416919 0.4932372520354784 --pwms 0.4101004931767698 0.642723289 0.4981727770042379 --pwmm 0.4101004931767698 0.642723289 0.4981727770042379 --pwss 0.4101004931767698 0.642723289 0.4981727770042379 \
--emsa -0.42745284 -0.25637422 0.0182113111 --emsn -1.15496282 -0.205819368 0.0287151085 \
--emma -13.7680371 -14.2357021 -14.1778536 --emmn -13.31637158 -13.3154217 -14.2030666 \
--essa -0.22232703 -0.24642711 -0.00560155865 --essn -0.36354996 0.155786444 0.039231938 \
--csv PIB-diisobutyl_ketone_227k.csv PIB-diisobutyl_ketone_285k.csv PIB-diisobutyl_ketone_6m.csv \
--img PIB-DIBK-cmaes_tot.png --dump PIB-DIBK-227k.dump PIB-DIBK-285k.dump PIB-DIBK-6m.dump

