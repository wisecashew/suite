#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time

if __name__=="__main__":

	start = time.time()

	settings = "pnipam_tip4p-ice-water.settings"
	topo     = "pnipam_tip4p-ice-water.r1.data"
	traj     = "coords.lammpstrj"
	dt       = 2

	my_settings = lmp_settings.Settings(settings)
	my_settings.read_settings()

	my_data     = lmp_data.Data(topo)
	my_data.read_data()

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	H_polymer = "type 8"
	H_water   = "type 15"
	O_polymer = "type 6"
	N_polymer = "type 7"
	O_water   = "type 14"

	# access atom types for hbonding
	print(f"Run hbond analysis...")
	hbonds_with_Opoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=O_polymer, donors_sel=O_water)
	hbonds_with_Npoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=N_polymer, donors_sel=O_water)
	hbonds_with_Owater = HBA(universe=u, hydrogens_sel=H_polymer, acceptors_sel=O_water,   donors_sel=N_polymer)

	# run hbond analysis with MDA
	hbonds_with_Opoly.run()
	hbonds_with_Npoly.run()
	hbonds_with_Owater.run()

	print(f"Done running.")
	print(f"Amidic N with Water H: \n{hbonds_with_Npoly.results.hbonds.shape}")
	print(f"Carbonyl O with Water H: \n{hbonds_with_Opoly.results.hbonds.shape}")
	print(f"Water O with Amide H: \n{hbonds_with_Owater.results.hbonds.shape}")
	
	print(f"Completed computation.")

	# fig = plt.figure(figsize=(2,2))
	# ax  = plt.axes()

	# ax.plot(times/1e+6, counts, lw=1, ls='--', color='steelblue', marker='o', markersize=2, mec='k')
	# fig.savefig(args.img, bbox_inches="tight", dpi=1200)

	# stop = time.time()
	# print(f"Time for computation is {stop-start} seconds.", flush=True)
