#!/bin/bash

set -e

for string in db1/chisc*.db; do
	echo $string
	chisc=$(echo "$string" | grep -oP '(?<=chisc_)[0-9.]+')
	chips=$(echo "$string" | grep -oP '(?<=chips_)[0-9.]+')
	chipc=$(echo "$string" | grep -oP '(?<=chipc_)[0-9.]+')
	vs=$(echo "$string" | grep -oP '(?<=vs_)[0-9.]+')
	vp=$(echo "$string" | grep -oP '(?<=vp_)[0-9.]+')
	vc=$(echo "$string" | grep -oP '(?<=vc_)[0-9.]+')
	python plot_db.py --db $string --img images1/chisc_${chisc}-chips_${chips}-chipc_${chipc}-vs_${vs}-vp_${vp}-vc_${vc}.png
done
