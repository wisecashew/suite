#!/bin/bash

set -e

chisc=1
chips=0
chipc=5
vs=1
vc=4
vp=3

python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/bcalculator.py \
--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp \
--crit-pkl crits/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl \
--mesh-pkl /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/mesh.pkl \
--final-binodal  chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.binodal.pkl --plot-crits --plot-binodal \
--island-stable-pkl stable_islands/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.stable_islands.pkl \
--island-unstable-pkl unstable_islands/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.unstable_islands.pkl 
