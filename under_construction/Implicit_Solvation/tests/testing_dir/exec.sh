#!/bin/bash
set -e 

rm stderr.txt || true

./MonteCarlo -f 1 -M 20000 -p positions.txt -u energydump.txt -t geom_and_esurf.txt -o coords.txt -v > stderr.txt 2>&1
echo "ran MC."
