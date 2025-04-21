#!/bin/awk -f

BEGIN {ignore=0}

$1 >= 48000000 && $1 <= 55990000 {ignore=1; next;}
$1 > 55990000                    {ignore=0;}
!ignore                          {print $0;}

END {}
