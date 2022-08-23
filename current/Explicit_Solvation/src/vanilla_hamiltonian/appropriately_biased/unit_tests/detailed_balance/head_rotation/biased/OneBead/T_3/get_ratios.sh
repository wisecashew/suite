#!/bin/bash 

T=`grep 'Temperature' o_to_n | grep -oE '\d+\.\d+|\d+'`
hits_on=`grep 'hits' o_to_n | grep -oE '\d+'`
hits_no=`grep 'hits' n_to_o | grep -oE '\d+'`
Eo=`grep 'Energy of system' o_to_n | grep -oE '\d+\.\d+|\d+'`
En=`grep 'Energy of system' n_to_o | grep -oE '\d+\.\d+|\d+'`
# b=$(( -1/$T * (En-Eo) ))

lhs=$(awk "BEGIN{print $hits_on/$hits_no }")
echo "Ratio of empirical transition probabilities = $lhs"

rhs=$(awk "BEGIN{print exp(1/$T * ($En-$Eo)) }")
echo "Ratio of configurational weights rho_n/rho_o = $rhs"

# echo "$T $hits_on $hits_no $Eo $En"

