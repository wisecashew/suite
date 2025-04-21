#!/bin/awk -f

BEGIN {ignore_config=0; nzero=0; dts=10000; current=0;}

/^ITEM: TIMESTEP$/                     {hold_line=$0; ignore_config=0; after_timestep=1; next}
/^0$/ && nzero>0                       {ignore_config=1;}
ignore_config                          {next;}
after_timestep && /^0$/ && nzero==0    {nzero+=1;     print hold_line; print current; after_timestep=0; next;}
after_timestep && /^[0-9]+$/ && !/^0$/ {current+=dts; print hold_line; print current; after_timestep=0; next;}
!ignore_config                         {print $0;}

END {}
