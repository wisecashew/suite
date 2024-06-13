#!/bin/bash

set -e
python nuanced_arr.py --label UCST --skips 500 --draw-bin \
--xlim 0 0.4 --ylim 280 340 \
--Tbot 200 200 200 --Ttop 400 400 400 --T 200 400 --Tcrit 275 350 \
--vs 0.94846692 1.07852990 1.08670771 --vm 150.9863077 999.823864 4324.01427 \
--pv 1.0 1.0 1.0 --pwms 0.255717 0.229323 0.244438 --pwmm 0 0 0 --pwss 0 0 0 \
--emsa -1.00333732 -1.05444885 -0.958422900 --emsn -0.9956613 -0.855031947 -1.04073626 \
--emma -16.99373763 -14.9983981 -15.0910125 --emmn -16.99373763 -14.9983981 -15.0910125 \
--essa 0 0 0 --essn 0 0 0 \
--csv PIB-diisobutyl_ketone_227k.csv PIB-diisobutyl_ketone_285k.csv PIB-diisobutyl_ketone_6m.csv \
--img PIB-DIBK-cmaes_tot.svg --dump PIB-DIBK-227k.dump PIB-DIBK-285k.dump PIB-DIBK-6m.dump

python nuanced_arr.py --label UCST --skips 500 --draw-bin --xlim 0 0.4 --ylim 280 340 --Tbot 200 200 200 --Ttop 400 400 400 --T 200 400 --Tcrit 275 350 --vs 0.94846692 1.07852990 1.08670771 --vm 150.9863077 999.823864 4324.01427 --pv 1.0 1.0 1.0 --pwms 0.255717 0.229323 0.244438 --pwmm 0 0 0 --pwss 0 0 0 --emsa -1.00333732 -1.05444885 -0.958422900 --emsn -0.9956613 -0.855031947 -1.04073626 --emma -16.99373763 -14.9983981 -15.0910125 --emmn -16.99373763 -14.9983981 -15.0910125 --essa 0 0 0 --essn 0 0 0 --csv PIB-diisobutyl_ketone_227k.csv PIB-diisobutyl_ketone_285k.csv PIB-diisobutyl_ketone_6m.csv --img PIB-DIBK-cmaes_tot.svg --dump PIB-DIBK-227k.dump PIB-DIBK-285k.dump PIB-DIBK-6m.dump
:<<'END'
python nuanced_arr.py --label UCST --skips 1000 --draw-bin \
--xlim 0 0.4 --ylim 280 340 \
--Tbot 200 --Ttop 400 --T 200 400 --Tcrit 275 350 \
--vs 1.08670771 --vm 4324.01427 \
--pv 1.0 --pwms 0.244438 --pwmm 0 --pwss 0 \
--emsa -0.958422900 --emsn -1.04073626 \
--emma -15.0910125 --emmn -15.0910125 \
--essa 0 --essn 0 \
--csv PIB-diisobutyl_ketone_6m.csv \
--img PIB-DIBK-cmaes_tot.png --dump PIB-DIBK-6m.dump
END
