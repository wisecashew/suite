#!/bin/bash

set -e

for f in BINODALS/*; do
	echo "In $f..."
	chisc=$(echo "$f" | grep -oP '(?<=chisc_)-?[0-9]+(\.[0-9]+)?')
	chips=$(echo "$f" | grep -oP '(?<=chips_)-?[0-9]+(\.[0-9]+)?')
	chipc=$(echo "$f" | grep -oP '(?<=chipc_)-?[0-9]+(\.[0-9]+)?')
	vs=$(echo "$f" | grep -oP '(?<=vs_)-?[0-9]+(\.[0-9]+)?')
	vc=$(echo "$f" | grep -oP '(?<=vc_)-?[0-9]+(\.[0-9]+)?')
	vp=$(echo "$f" | grep -oP '(?<=vp_)-?[0-9]+(\.[0-9]+)?')
	bi="chips_${chips}-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl"
	ci="chips_${chips}-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl"
	img="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.png"
	python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/prober.py \
	--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp --mesh-pkl ../mesh.pkl \
	--binodal-pkl BINODALS/$bi --crit-pkl CRITS/$ci --img IMAGES/$img --plot-crits 
done
