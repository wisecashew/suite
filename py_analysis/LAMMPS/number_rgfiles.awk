#!/bin/awk -f

BEGIN {current=0; dt=5000;}

1{$1=current; print; current+=dt}

END {}
