#!/bin/bash

# ./flogotts -f 1 -M 100 -l 1000 -p 8mer.txt -t geom_and_esurf.txt -u energy_.mc -o coords_.mc -e orientations_.mc -s stats_.mc -H solvation_shell_.mc -L latdump_.mc -R latdump_.mc -r -v

# ./flogotts -f 1 -M 30000 -l 1000 -p 8mer.txt -t geom_and_esurf.txt -u energy.mc -o coords.mc -e orientations.mc -s stats.mc -H solvation_shell.mc -L latdump.mc --isotropic

# ./flogotts -f 1 -M 1000 -l 1000 -p 8mer.txt -t fhp.txt -u energy.mc -o coords.mc -e orientations.mc -s stats.mc -L latdump.mc
./flogotts -f 1 -M 1000 -l 1000 -p 8mer.txt -t fhp.txt -u energy.mc -o coords.mc -e orientations.mc -s stats.mc -L latdump.mc -R latdump.mc -r
