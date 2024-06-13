#!/bin/bash

set -e

for f in binodals/*; do
	echo "In $f..."
	chisc=$(echo "$f" | grep -oP '(?<=chisc_)-?[0-9]+(\.[0-9]+)?')
	chips=$(echo "$f" | grep -oP '(?<=chips_)-?[0-9]+(\.[0-9]+)?')
	chipc=$(echo "$f" | grep -oP '(?<=chipc_)-?[0-9]+(\.[0-9]+)?')
	vs=$(echo "$f" | grep -oP '(?<=vs_)-?[0-9]+(\.[0-9]+)?')
	vc=$(echo "$f" | grep -oP '(?<=vc_)-?[0-9]+(\.[0-9]+)?')
	vp=$(echo "$f" | grep -oP '(?<=vp_)-?[0-9]+(\.[0-9]+)?')
	bi="chips_${chips}-chipc_${chipc}-chisc_${chisc}-vs_${vs}-vc_${vc}-vp_${vp}.binodals.pkl"
	# chisc_=$(echo "$chisc" | sed 's/\.0$//')
	# chips_=$(echo "$chips" | sed 's/\.0$//')
	# chipc_=$(echo "$chipc" | sed 's/\.0$//')
	vs_=$(echo "$vs" | sed 's/\.0$//')
	vc_=$(echo "$vc" | sed 's/\.0$//')
	vp_=$(echo "$vp" | sed 's/\.0$//')
	ci="chips_${chips}-chipc_${chipc}-chisc_${chisc}-vs_${vs_}-vc_${vc_}-vp_${vp_}.crits.pkl"
	img="chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.png"
	echo "$ci"
	python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/prober.py \
	--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp --mesh-pkl ../mesh.pkl \
	--binodal-pkl binodals/$bi --crit-pkl crits/$ci --img images/$img --plot-crits 
done
