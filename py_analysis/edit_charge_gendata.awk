#!/bin/awk -f
BEGIN {atoms=0; Onum=10; Hnum=11;}

!atoms && /OW/ {Onum=$1}
!atoms && /HW/ {Hnum=$1}

/^Atoms$/ {atoms=1}
atoms && $3==Onum {$4="-0.8476"}
atoms && $3==Hnum {$4="0.4238"}
/^Bonds$/ {atoms=0}
1 {print}

END {}
