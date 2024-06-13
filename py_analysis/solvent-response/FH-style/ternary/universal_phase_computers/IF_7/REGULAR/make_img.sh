#!/bin/bash

set -e

for f in more_binodal/*; do
	echo "In $f..."
	chisc=$(echo "$f" | grep -oP 'chisc_\K\d+')
	chips=$(echo "$f" | grep -oP 'chips_\K\d+')
	chipc=$(echo "$f" | grep -oP 'chipc_\K\d+')
	vs=$(echo "$f" | grep -oP 'vs_\K\d+')
	vc=$(echo "$f" | grep -oP 'vc_\K\d+')
	vp=$(echo "$f" | grep -oP 'vp_\K\d+')
	bi="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl"
	ci="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl"
	img="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.png"
	python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/binodal_parser.py \
	--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp \
	--final-binodal more_binodals/$bi --crit-pkl crits/$ci --img images/$img --plot-crits 
done
