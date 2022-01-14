#!/bin/bash 
set -e 

for T in 0.1 0.2 0.25 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 
do 
	cd Tequal${T}
	cp heat_cap.txt ../heat_capacity/heat_cap_$T.txt
	cp energy_freq.png ../energy_freq/energy_freq_$T.png
	cp rg_kT$T.png ../rg_hist/.
	cp stat_kt$T.txt ../rg_hist/. 
	cd .. 
done 

