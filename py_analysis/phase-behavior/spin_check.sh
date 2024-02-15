#!/bin/bash

set -e

python spin_plotter.py -pv 1 \
-emma -2500 \
-emmn -2000 \
-emsa -1900 \
-emsn -200 \
-essa -1000 \
-essn 0 \
-pwms 0.2 \
-pwmm 0.25 \
-pwss 0.4 \
-colors "coral" \
--png-name random_check
