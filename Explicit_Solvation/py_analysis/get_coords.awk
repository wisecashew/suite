#!/bin/awk -f

BEGIN{ 
    print "START POLYMER 0"
    regex="Dumping coordinates at step "step_var
    FS="|"
}

$0 ~ regex { f=1; }
print_flag && /END/ { exit; }
print_flag { print $1 " " $2 " " $3 }
f && /START/ { print_flag=1; f=0; }

END {
    print "END POLYMER 0" 
}
