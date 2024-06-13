#!/bin/bash

set -e


for file in chisc*..out; do
	# if grep -q "Detected 1 oom_kill" $file; then 
	if grep -q "srun: error" $file; then 
		echo "In $file..."
		chisc=$(echo "$file" | grep -oP '(?<=chisc_)-?\d+(\.\d+)?' | sed 's/chisc_//')
		chips=$(echo "$file" | grep -oP '(?<=chips_)-?\d+(\.\d+)?' | sed 's/chips_//')
		chipc=$(echo "$file" | grep -oP '(?<=chipc_)-?\d+(\.\d+)?' | sed 's/chipc_//')
		vs=$(echo "$file" | grep -oP '(?<=vs_)-?\d+(\.\d+)?' | sed 's/vs_//')
		vp=$(echo "$file" | grep -oP '(?<=vp_)-?\d+(\.\d+)?' | sed 's/vp_//')
		vc=$(echo "$file" | grep -oP '(?<=vc_)-?\d+(\.\d+)?' | sed 's/vc_//')
		echo "f = $file: chisc = $chisc, chips = $chips, chipc = $chipc, vs = $vs, vp = $vp, vc = $vc"
	fi
done
