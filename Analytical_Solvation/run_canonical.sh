#!/bin/bash

set -e

# python canonical_partition-v2.py -o FHP-flat --EMMA -2     --EMMN -2     --EMSA -1 --EMSN -1 --color "#369DE8"
# python canonical_partition-v2.py -o FHP-R0 --EMMA -1     --EMMN -1     --EMSA -1 --EMSN -1 --color "#369DE8"
# python canonical_partition-v2.py -o FHP-R1 --EMMA -3     --EMMN -1     --EMSA -1 --EMSN -1 --color "#1FB967"
# python canonical_partition-v2.py -o FHP-R2 --EMMA -1.5   --EMMN -1.5   --EMSA -1 --EMSN 0  --color "#B9B41F"
python canonical_partition-v2.py -o FHP-R3-p1 --EMMA -2.001  --EMMN -2.001  --EMSA -1 --EMSN 0  --color "#B91F72"
python canonical_partition-v2.py -o FHP-R3-p2 --EMMA -2.0001 --EMMN -2.0001 --EMSA -1 --EMSN 0  --color "#B91F72"
