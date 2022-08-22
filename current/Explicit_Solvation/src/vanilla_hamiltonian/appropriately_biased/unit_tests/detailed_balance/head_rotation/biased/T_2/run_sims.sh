#!/bin/bash 

echo "running from initial to final..."
./hr_o_n_boltz -f 10000 -M 1'000'000 -p initial.txt -t geom_and_esurf.txt -u energy -L lattice -e orientation -s stats -o coords >> o_to_n 

echo "running from final to initial..."
./hr_n_o_boltz -f 10000 -M 1'000'000 -p final.txt -t geom_and_esurf.txt -u energy -L lattice -e orientation -s stats -o coords >> n_to_o
