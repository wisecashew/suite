#!/bin/bash

set -e

./../McLattE.exe fhp -f 1 -M 1000 -l 100 -t geom_and_esurf.txt -o coords.mc \
-p equi32mer_1.txt -u energy.mc -s stats.mc -e orientation.mc \
-L latdump.mc -H solvation.mc
