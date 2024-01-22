#!/bin/bash

set -e

# Loop through files with ".unstable_points." in the name
for file in UNSTABLE_ISLANDS/chips_*chipc_*chisc_*vs*vc*vp*.unstable_islands.*; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Add your processing logic here	
		chips=$(echo "$file" | grep -oP 'chips_-?[\d.]+' | sed 's/chips_//')
		chipc=$(echo "$file" | grep -oP 'chipc_-?[\d.]+' | sed 's/chipc_//')
		chisc=$(echo "$file" | grep -oP 'chisc_-?[\d.]+' | sed 's/chisc_//')
		vs=$(echo "$file" | grep -oP 'vs_-?[\d.]+' | sed 's/vs_//')
		vc=$(echo "$file" | grep -oP 'vc_-?[\d.]+' | sed 's/vc_//')
		vp=$(echo "$file" | grep -oP 'vp_-?\d+(\.\d+)?' | sed 's/vp_//')
		echo "Processing file: $file"
		echo "chips = $chips, chipc = $chipc, chisc = $chisc, vs = $vs, vc = $vc, vp = $vp"

		# srun --ntasks=1 --nodes=1 --cpus-per-task=1 python universal.py --chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp --island-stable-pkl STABLE_ISLANDS/chips_$chips-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.stable_islands.pkl --island-unstable-pkl UNSTABLE_ISLANDS/chips_$chips-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.unstable_islands.pkl --crit-pkl CRITS/chips_$chips-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl --binodal-pkl BINODALS/chips_$chips-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl --search-density 1000 --final-binodal BINODALS/chips_$chips-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl --plot-binodals --plot-edges --plot-crits > big_search.chips_${chips}-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.out 2>&1 &
        # Add more commands as needed
    else
        echo "File not found: $file"
    fi
done

