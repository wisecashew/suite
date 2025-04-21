#!/bin/awk -f

BEGIN {dts=10000; current=0;}

/^ITEM: TIMESTEP$/           {hold_line=$0; after_timestep=1; next;}
after_timestep && /^[0-9]+$/ {print hold_line; print current; current+=dts; after_timestep=0; next;}
!after_timestep              {print; next;}

END {}
