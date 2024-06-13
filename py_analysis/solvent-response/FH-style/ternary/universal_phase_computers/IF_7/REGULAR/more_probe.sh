#!/bin/bash

set -e

chisc=0
chips=1
chipc=4
vs=1
vc=1
vp=4

python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/prober.py \
--chisc $chisc --chips $chips --chipc $chipc -vs $vs -vc $vc -vp $vp \
--crit-pkl crits/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl \
--binodal-pkl chisc_${chisc}_chips_${chipc}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.binodal.pkl --plot-crits --plot-binodal 
