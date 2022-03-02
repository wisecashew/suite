#!/bin/bash 

./SimpMC -f 1 -M 10000 -p positions.txt -u energydump.txt -t geom_and_esurf.txt -o coords.txt > stderr.txt 2>&1
echo "executed simulation..."

grep '[0-9]+\.5' energydump.txt 
echo "grepped through energydump.txt" 
grep "Time taken " stderr.txt 
grep "Number of moves accepted" stderr.txt 
