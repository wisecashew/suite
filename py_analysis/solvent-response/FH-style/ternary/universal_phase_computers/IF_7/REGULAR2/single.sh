#!/bin/bash

set -e

chisc=0
chips=1
chipc=3
vs=5
vc=4
vp=1


python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/bcalculator.py --chisc ${chisc} --chips ${chips} --chipc ${chipc} -vs ${vs} -vc ${vc} -vp ${vp} --final-binodal binodals/chisc_${chisc}_chips_${chips}_chipc_${chipc}_vs_${vs}_vc_${vc}_vp_${vp}.pkl --island-stable-pkl stable_islands/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.stable_islands.pkl --island-unstable-pkl unstable_islands/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.unstable_islands.pkl --mesh-pkl ../../mesh.pkl --crit-pkl crits/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl --plot-binodals

# python /scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/solvent-response/FH-style/ternary/universal_phase_computers/prober.py --chisc ${chisc} --chips ${chips} --chipc ${chipc} -vs ${vs} -vc ${vc} -vp ${vp} --binodal-pkl binodals/chisc_${chisc}_chips_${chips}_chipc_${chipc}_vs_${vs}_vc_${vc}_vp_${vp}.pkl --mesh-pkl ../../mesh.pkl --crit-pkl crits/chisc_${chisc}_chips_${chips}_chipc_${chipc}-vs_${vs}-vc_${vc}-vp_${vp}.crits.pkl --plot-binodals --img my_bigmem
