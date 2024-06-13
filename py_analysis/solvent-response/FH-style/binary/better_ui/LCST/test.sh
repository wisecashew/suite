#!/bin/bash

set -e

python ../regularized-fits.py --label UCST --T 400 500 \
--T-UCST 400 --Tbot 400 --Ttop 500 \
--csv PS-ethylacetate_100k.csv --sig 0.01 --L1 0.00 --L2 0.00 \
--skip 5000 --img test_100k --inp PS-EA-100k_.in

:<<'END'
python ../regularized-fits.py --label LOOP --T 400 500 \
--T-loop 350 --Tbot 400 --Ttop 500 \
--csv PS-ethylacetate_223k.csv --sig 0.01 --L1 0.00000000001 --L2 0.00000000001  \
--skip 500 --img test_223k --inp PS-EA-223k_.in

python ../regularized-fits.py --label LOOP --T 400 500 \
--T-loop 350 --Tbot 400 --Ttop 500 \
--csv PS-ethylacetate_600k.csv --sig 0.01  --L1 0.00000000001 --L2 0.00000000001 \
--skip 500 --img test_600k --inp PS-EA-600k_.in
END
