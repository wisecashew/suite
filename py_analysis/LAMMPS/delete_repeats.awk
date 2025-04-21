#!/bin/awk -f

BEGIN {nzero=0;}

/^ITEM: TIMESTEP$/                     {hold_line=$0; next;}
/^0$/                                  {nzero+=1}
nzero == 8                             {exit;}
/^0$/ && nzero < 8                     {print hold_line; print $0; next;}
1                                      {print; next;}

END {}
