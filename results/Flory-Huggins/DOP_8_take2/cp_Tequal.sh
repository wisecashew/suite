#!/bin/bash 

for T in 0.1 0.2 0.25 0.3 0.4 0.45 
do
	cp -r Tequal0.5 Tequal${T}
	cd Tequal${T} 
	sed "s/kT = 0.5/kT = ${T}/g" geom_and_esurf.txt > tmp 
	mv tmp geom_and_esurf.txt 
	cd .. 
done

