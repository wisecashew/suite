#!/bin/bash 

for T in 0.1 0.2 0.25 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 
do
	cp py_analysis/{PED.py,PRG.py} Tequal$T/.
done
echo "sent!"
	
