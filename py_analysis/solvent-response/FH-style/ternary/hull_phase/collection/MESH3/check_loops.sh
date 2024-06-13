#!/bin/bash

set -e

count=0

for vs in 2 2.5 3; do
	for vc in 2 2.5 3; do
		for vp in 2 2.5 3; do
			for chisc in 1 1.5 2; do
				for chips in 1 1.5 2; do
					for chipc in 1 1.5 2; do
						if ([ "$vs" != "$vc" ] || [ "$vs" != "$vp" ]) && ([ "$chisc" != "$chips" ] || [ "$chipc" != "$chisc" ]); then
							count=$(($count+1))
							echo "vs = $vs, vc = $vc, vp = $vp, chisc = $chisc, chips = $chips, chipc = $chipc"
						fi
					done
				done
			done
		done
	done
done

echo "Number of parameters = $count."
