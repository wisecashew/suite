#!/bin/bash

set -e


python ../regularized-fits.py --T 400 600 --T-loop 475 \
--inp PEG-water.in  --sig 0.1 \
--img PEG-water --csv PEG-water.csv \
--Tbot 400 --Ttop 550 --label LOOP --skip 500 --L1 0.00001 --L2 0.00001 

python ../regularized-fits.py --T 300 410 --T-loop 350 \
--inp 2BE-water.in --sig 0.1 \
--img 2BE-water --csv 2-butoxyethanol-water.csv \
--Tbot 300 --Ttop 400 --label LOOP --skip 500 --L1 0.00001 --L2 0.00001 
