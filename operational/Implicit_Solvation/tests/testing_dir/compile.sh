#!/bin/bash

./MonteCarlo -f 1 -M 20000 -p positions.txt -u energydump.txt -t geom_and_esurf.txt -o coords.txt -v > stderr.txt 2>&1

echo "ran MC. "
grep '[0-9]\.5\|keys is 7' stderr.txt
grep '[0-9]\.5' energydump.txt 
echo "ran grep."
