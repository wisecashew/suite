#!/bin/bash 

for dir in *
do
	if [ -d "$dir" ]; then 
		cd $dir 
		rm energy_freq.png stat_kt0.5.txt rg_kT0.5.png MonteCarlo
		cd .. 
	fi
done  
