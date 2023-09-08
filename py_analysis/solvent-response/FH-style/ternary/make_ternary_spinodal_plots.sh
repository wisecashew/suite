#!/bin/bash

set -e

for chisc in -11; do
	python cspin.py --chisc ${chisc} --chips -1 --chipc -1 -N 32 --ternary &
done
# wait

for chipc in -{1..10}; do
	python cspin.py --chisc 0 --chips -1 --chipc ${chipc} -N 32 --ternary &
done
wait
