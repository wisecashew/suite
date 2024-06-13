#!/bin/bash

set -e


python constricted-cmaes-with_labels.py --label LOOP --T 400 550 --T-loop 450 --skip 1000 --vm 141.7021338078609 \
--vs 7.691265982329057 --xfrmpv 23.10019043182071 --xfrmpw -1.2368622351827026 --emm -153.12454913378573 --emsa -127.02258357518589 --emsn -56.94316370374701 --ess 0 \
--csv PEG-water.csv --img peg-water-cons --optimize --maxiter 100
