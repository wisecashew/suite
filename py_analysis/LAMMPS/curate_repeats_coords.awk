#!/bin/awk -f

BEGIN {ignore_config=0;}

/^ITEM: TIMESTEP$/                                 {hold_line=$0; ignore_config=0; after_timestep=1; next}
$1 >= 48000000 && $1 <= 55980000 && after_timestep {ignore_config=1; next;}
ignore_config                                      {next;}
!ignore_config && after_timestep                   {print hold_line; print $0; after_timestep=0; next;}
!ignore_config                                     {print $0; next;}

END {}
