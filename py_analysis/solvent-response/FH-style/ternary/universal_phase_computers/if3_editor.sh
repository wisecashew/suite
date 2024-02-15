#!/bin/bash

set -e

cd IF_3
cd stable_islands
for i in 2.{4..9} 3.0; do
	for j in 2.{4..9} 3.0; do
		for k in 2.{4..9} 3.0; do
			cp chips_${i}-chipc_${j}-chisc_${k}-vs_1-vc_1-vp_1.stable_islands.pkl ../stable_islands_copy/chips_${j}-chipc_${k}-chisc_${i}-vs_1-vc_1-vp_1.stable_islands.pkl 
		done
	done
done

cd ../unstable_islands

for i in 2.{4..9} 3.0; do
	for j in 2.{4..9} 3.0; do
		for k in 2.{4..9} 3.0; do
			cp chips_${i}-chipc_${j}-chisc_${k}-vs_1-vc_1-vp_1.unstable_islands.pkl ../unstable_islands_copy/chips_${j}-chipc_${k}-chisc_${i}-vs_1-vc_1-vp_1.unstable_islands.pkl 
		done
	done
done



