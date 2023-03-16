#!/bin/bash 

set -e 

for sim in biased unbiased; do 
    for T in 1 2 3; do
        cp geom_and_esurf.txt $sim/T_$T/.
        cp 8mer.txt $sim/T_$T/. 
    done
done 
 
