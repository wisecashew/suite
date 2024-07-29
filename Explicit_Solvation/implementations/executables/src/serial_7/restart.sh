#!/bin/bash

# ./flogotts -f 1 -M 100 -l 1000 -p 8mer.txt -t geom_and_esurf.txt -u energy_.mc -o coords_.mc -e orientations_.mc -s stats_.mc -H solvation_shell_.mc -L latdump_.mc -R latdump_.mc -r -v

./flogotts -f 1 -M 10000 -l 1000 -p 8mer.txt -t geom_and_esurf.txt -u energy_.mc -o coords.mc -e orientations_.mc -s stats_.mc -H solvation_shell_.mc -L latdump.mc -R latdump.mc -r -v

# McLattE -f 1 -M 5112 -l 1000 -p 8mer.txt -t geom_and_esurf.txt -u energy.mc -o coords.mc -e orientations.mc -s stats.mc -L latdump.mc
