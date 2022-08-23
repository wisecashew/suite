#!/bin/bash 

echo "running from initial to final..."
./hr_o_n_boltz -f 100000 -M 10'000'000 -p initial.txt -t geom_and_esurf.txt -u energy1 -L lattice1 -e orientation1 -s stats1 -o coords1 > o_to_n & 

echo "running from final to initial..."
./hr_n_o_boltz -f 100000 -M 10'000'000 -p final.txt -t geom_and_esurf.txt -u energy2 -L lattice2 -e orientation2 -s stats2 -o coords2 > n_to_o & 

wait 
