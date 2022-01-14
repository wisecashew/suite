#!/bin/bash 

for T in 0.1 0.2 0.25 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 
do
	cd Tequal${T}
	echo "entered T = $T."
	python PRG.py -i coords.txt -e 10 -T $T
	echo "calculated Rg..." 
	python PED.py -i energydump.txt -eh energy_freq.png -c heat_cap.txt -T $T
	echo "ran python in file Tequal${T}..." 
	cd ..
done 

cd heat_capacity 
python Cv.py
cd ..

echo "Executed all python analysis."

