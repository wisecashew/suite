#!/bin/bash

set -e

python ../regularized-fits.py --label UCST --T 400 500 \
--T-UCST 400 --Tbot 400 --Ttop 500 \
--csv PS-ethylacetate_100k.csv --sig 0.01 --L1 0.00 --L2 0.00 \
--skip 5000 --img PS-EA-100k --maxiter 10 --optimize --inp PS-EA-100k_.in --oup PS-EA-100k_.in
