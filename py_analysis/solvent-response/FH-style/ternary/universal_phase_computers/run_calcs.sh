#!/bin/bash

set -e

for c in 2.4 2.5 2.6 2.7 2.8 3.0; do

	python compute.py --chisc $c --chips $c --chipc $c -vs 1 -vc 1 -vp 1 --plot-crits 

done
