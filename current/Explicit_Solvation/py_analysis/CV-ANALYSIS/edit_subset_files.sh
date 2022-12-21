#!/bin/bash 

fnames=`ls *subset.py`

for f in $fnames; do
	echo "In $f..."
	# sed 's/.*U_list = .*/    U_list = ["U8", "U9", "U10"]/' $f > sedput 
	sed 's/-2000/-1500/g' $f > sedput
	mv sedput ${f}
done

chmod +x *.py
