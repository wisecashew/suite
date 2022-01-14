#!/bin/bash 

for T in 1 2 3 4 5
do
	cd Tequal$T
	rm MonteCarlo
	cd ..
	cp MonteCarlo Tequal$T/.
done
