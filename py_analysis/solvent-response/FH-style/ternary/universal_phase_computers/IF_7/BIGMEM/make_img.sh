#!/bin/bash

set -e

for f in binodals/*; do

	if echo $f | grep -qE "\.pkl$"; then
		echo "In $f..."
		chisc=$(echo "$f" | grep -oP '(?<=chisc_)\d+\.\d+')
		chips=$(echo "$f" | grep -oP '(?<=chips_)\d+\.\d+')
		chipc=$(echo "$f" | grep -oP '(?<=chipc_)\d+\.\d+')
		vs=$(echo "$f" | grep -oP '(?<=vs_)\d+\.\d+')
		vc=$(echo "$f" | grep -oP '(?<=vc_)\d+\.\d+')
		vp=$(echo "$f" | grep -oP '(?<=vp_)\d+\.\d+')
		bi="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl"
		chisc_=$(echo "$chisc" | sed 's/\.0$//')
		chips_=$(echo "$chips" | sed 's/\.0$//')
		chipc_=$(echo "$chipc" | sed 's/\.0$//')
		vs_=$(echo "$vs" | sed 's/\.0$//')
		vc_=$(echo "$vc" | sed 's/\.0$//')
		vp_=$(echo "$vp" | sed 's/\.0$//')
		ci="chisc_${chisc_}_chips_${chips_}_chipc_${chipc_}-vs_${vs_}-vc_${vc_}-vp_${vp_}.crits.pkl"
		img="chisc_${chisc_}_chips_${chips_}_chipc_${chipc_}-vs_${vs_}-vc_${vc_}-vp_${vp_}.png"
		python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/prober.py \
		--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp \
		--binodal-pkl binodals/$bi --crit-pkl crits/$ci --img images/$img --plot-crits 
	else
		echo "$f -- moving on."
	fi
done
