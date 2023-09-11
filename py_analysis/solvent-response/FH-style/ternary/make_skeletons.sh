#!/bin/bash

set -e


echo "Inside block 1..."
for vp in {2..10..2}; do
	for chisc in -{8..12}; do
		echo "vp = $vp, chisc = $chisc..."
		python pskelbin_uni2.py --chisc $chisc --chips -1 --chipc -1 -vs 1 -vc 1 -vp $vp --mesh 200 --no-rtw > chisc_$chisc-vp_$vp.out 2>&1
	done
done




