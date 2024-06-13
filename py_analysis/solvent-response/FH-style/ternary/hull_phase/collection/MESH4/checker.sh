#!/bin/bash

set -e

for file in chisc*.out; do
	# echo "In $file..."
	chisc=$(echo "$file" | grep -oP '(?<=chisc_)-?\d+(\.\d+)?' | sed 's/chisc_//')
	chips=$(echo "$file" | grep -oP '(?<=chips_)-?\d+(\.\d+)?' | sed 's/chips_//')
	chipc=$(echo "$file" | grep -oP '(?<=chipc_)-?\d+(\.\d+)?' | sed 's/chipc_//')
	vs=$(echo "$file" | grep -oP '(?<=vs_)-?\d+(\.\d+)?' | sed 's/vs_//')
	vp=$(echo "$file" | grep -oP '(?<=vp_)-?\d+(\.\d+)?' | sed 's/vp_//')
	vc=$(echo "$file" | grep -oP '(?<=vc_)-?\d+(\.\d+)?' | sed 's/vc_//')
	if [ -f chisc_${chisc}-chips_${chips}-chipc_${chipc}-vs_$vs-vp_$vp-vc_$vc.hull.pkl ]; then
		# echo "$file and conjugate exists."
		y=1
	else
		echo "$file and conjugate DO NOT exist."
	fi
done
