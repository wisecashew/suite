#!/bin/bash

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

j=0
for file in collection/MESH2/*.pkl; do
	echo "$j: $file"
	chipc=$(echo "$file" | sed -n 's/.*chipc_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	chips=$(echo "$file" | sed -n 's/.*chips_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	chisc=$(echo "$file" | sed -n 's/.*chisc_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	vs=$(echo "$file" | sed -n 's/.*vs_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	vc=$(echo "$file" | sed -n 's/.*vc_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	vp=$(echo "$file" | sed -n 's/.*vp_\(-\?[0-9]*\.*[0-9]*\).*/\1/p')
	echo "$j: chisc = $chisc, chips = $chips, chipc = $chipc, vs = $vs, vc = $vc, vp = $vp"
	j=$(($j + 1))
done

