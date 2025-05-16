#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt
import field 

if __name__=="__main__":
	
	# Define the input parameters
	inputs = {
		"v_pp"     : -0.13778511,
		"v_pa"     : 0.08389742,
		"v_pb"     : 0.0,
		"v_aa"     : 0.0,
		"v_bb"     : 0.0,
		"v_ab"     : -0.081982783,
		"a_p"      : 1.0,
		"a_a"      : 0.1,
		"a_b"      : 0.1,
		"rho0"     : 1.1,
		"deg_poly" : 50,
	}

	# get the field equations for the homopolymer in dumbell solvent
	fld = field.Field(inputs)
	fld.print()

	# set up some figures for the diagram
	fig, ax = plt.subplots(1, 1, figsize=(3, 3))

	# get the mass fraction vector
	phi_range = np.hstack([[0], np.logspace(-6, np.log10(0.99), 100)])
	
	for p in phi_range:
		rho_p = fld.rho0 * p
		rho_s = fld.rho0 * (1-p)
		roots = fld.stab_roots(rho_p, rho_s)

		if len(roots) > 0:
			print(f"phi = {p}; T = {[1/r for r in roots]}", flush=True)
		
		if len(roots) == 2:
			for root in roots:
				ax.scatter([p], [1/root], c="steelblue", marker='o', ec='k', s=6)
		
		if len(roots) == 1:
			for root in roots:
				ax.scatter([p], [1/root], c="darkred", marker='o', ec='k', s=6)


	# set the axis limits
	ax.set_xlim(0, 1)
	# ax.set_xscale("log")
	# ax.set_yscale("log")
	fig.savefig("phase_diagram.png", dpi=1200)
