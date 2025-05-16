#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt
import field 

if __name__=="__main__":
	
	# Define the input parameters
	inputs = {
		"v_pp"     : -4.11215026,
		"v_pa"     : 0.62381171,
		"v_pb"     : 0.0,
		"v_aa"     : -0.13914001,
		"v_bb"     : -0.13914001,
		"v_ab"     :   0,
		"a_p"      : 1.0,
		"a_a"      : 0.1,
		"a_b"      : 0.1,
		"rho0"     : 5.1,
		"deg_poly" : 50,
	}

	# get the field equations for the homopolymer in dumbell solvent
	fld = field.Field(inputs)
	fld.print()

	# set up some figures for the diagram
	fig, ax = plt.subplots(1, 1, figsize=(3, 3))

	# set the densities
	rho_p = 0.1
	rho_s = 5

	# get the zero-temperatures
	roots = fld.stab_roots(rho_p, rho_s)
	print(f"thermodynamic beta roots = {[1/root for root in roots]}.", flush=True)

	# get the stability over temperature
	stab = fld.stability(rho_p, rho_s)

	print(f"This system is unstable at certain points: {np.any(stab < 0)}.", flush=True)
	ax.scatter(fld.T, stab, c="steelblue", marker='o', ec='k', s=10)
	ax.plot(fld.T, stab, c="steelblue", lw=0.5, ls='--')

	# set the axis limits
	ax.axhline(0, color='k', lw=0.5)
	ax.set_ylim(-0.5, 1)
	ax.set_xlim(0.01, 1e+2)
	ax.set_xscale("log")
	fig.savefig("stab.png", dpi=1200)
