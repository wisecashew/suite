#!/bin/bash
set -e 

shopt -s expand_aliases
source ~/.bash_profile

g++ $1 misc.cpp classes.cpp -o MonteCarlo

rm stderr.txt || true

./MonteCarlo -f 1 -N 1 -p positions.txt -t geom_and_esurf.txt > stderr.txt 2>&1
echo "ran MC."
