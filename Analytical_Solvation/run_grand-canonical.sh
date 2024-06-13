#!/bin/bash

set -e

python grand_canonical_partition-v3_.py -o PREFADS-v3 --EMM -1 --EMS -0.541666 --EMC -0.916666 --ESC 0 --color "steelblue" "coral" "forestgreen" "violet" "darkred" --markers "o" "^" "s" "*" "d" --T 0.1 0.2 0.5 1.0

python grand_canonical_partition-v3_.py -o PREFADS-v3 --EMM -1 --EMS -0.541666 --EMC -0.916666 --ESC 0 --color "steelblue" "coral" "forestgreen" "violet" "darkred" --markers "o" "^" "s" "*" "d" --T 0.1 0.2 0.5 1.0

python grand_canonical_partition-v3_.py -o SMIX-v3 --EMM -1 --EMS -0.541666 --EMC -1.541666 --ESC 0 --color "steelblue" "coral" "forestgreen" --markers "o" "^" "s" --T 0.1 1.0 10.0
