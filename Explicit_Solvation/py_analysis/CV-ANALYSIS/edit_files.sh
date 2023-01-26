#!/bin/bash 

fnames=`ls *T.py`

# U_list = aux.dir2U ( os.listdir (".") )
for f in $fnames; do
	echo "In $f..."
	# sed 's/.*U_list = .*/    U_list = aux.dir2U ( os.listdir(".") )/' $f > sedput 
	# sed 's/-2000/-1500/g' $f > sedput
	sed 's/.*divnorm = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ).*/        divnorm = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 )/' $f > sedput
	mv sedput ${f}
done

chmod +x *.py
