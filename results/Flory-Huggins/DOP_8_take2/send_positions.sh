#!/bin/bash 
set -e 

for T in 0.2 0.25 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 
do 
	
	cp Tequal0.1/positions_8.txt Tequal${T}/.
	
done 

