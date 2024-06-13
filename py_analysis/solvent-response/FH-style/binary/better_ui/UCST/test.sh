#!/bin/bash

set -e

python ../regularized-fits.py --label UCST --T 250 350 \
--T-UCST 250.0 --Tbot 250 --Ttop 350 \
--csv PIB-diisobutyl_ketone_227k.csv \
--skip 500 --img test_227k \
--inp PIB-DIBK-227k_.in --L1 0.01 --L2 0.01 > PIB-DIBK-227k.out 2>&1 

python ../regularized-fits.py --label UCST --T 250 350 \
--T-UCST 250.0 --Tbot 250 --Ttop 350 \
--csv PIB-diisobutyl_ketone_285k.csv \
--skip 500 --img test_285k \
--inp PIB-DIBK-285k_.in --oup PIB-DIBK-285k_.in --L1 0.01 --L2 0.01 > PIB-DIBK-285k.out 2>&1 

python ../regularized-fits.py --label UCST --T 250 350 \
--T-UCST 250.0 --Tbot 250 --Ttop 350 \
--csv PIB-diisobutyl_ketone_6m.csv \
--skip 500 --img test_6m \
--inp PIB-DIBK-6m_.in --L1 0.01 --L2 0.01 > PIB-DIBK-6m.out 2>&1 


