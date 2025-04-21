#!/home/satyend/.conda/envs/phase/bin/python

import matplotlib.pyplot as plt
import numpy as np
import time
import argparse
import pickle
import MDAnalysis.transformations as xforms
from mpl_toolkits.mplot3d import Axes3D
import latt_lmp as ll

parser = argparse.ArgumentParser(description="Compute lattice contact information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--Lpkl", dest="lpkl", type=str, action="store", help="enter address of processed lattice.", required=True)
parser.add_argument("--npkl", dest="npkl", type=str, action="store", help="enter address of new processed lattice.", required=True)
args = parser.parse_args()

if __name__=="__main__":

	# start the clock
	start = time.time()

	# extract the lattice object
	with open(args.lpkl, "rb") as latt_file:
		LL = pickle.load(latt_file)

	# calculate contacts
	for ts in LL.lattice:
		print(f"Calculating contacts @ timestep = {ts}. Calculating Nmm...", end=' ', flush=True)
		LL.get_Nmm(ts)
		print(f"done!", end=' ', flush=True)
		print(f"Calculating Nmsa...", end=' ', flush=True)
		LL.get_Nmsa(ts)
		print(f"done!", flush=True)

	# set up the new lattice object
	with open(args.npkl, "wb") as nfile:
		pickle.dump(LL.lattice, nfile)

	# nfile = open(args.npkl, "wb")
	# pickle.dump(LL, nfile)
	# nfile.close()

	"""
	# set up the figure
	fig = plt.figure(figsize=(3,3))
	ax  = fig.add_subplot(111, projection="3d")

	# define some timepoints
	timestep = 1000
	polymer  = LL.lattice[timestep]["polymer"]
	solvents = LL.lattice[timestep]["solvent"]
	msize = 1
	ssize = 1
	for monomer in polymer:
		ll.plot_cube(ax, monomer, msize, "coral")

	for solvent in solvents:
		ll.plot_cube(ax, solvent, ssize, "steelblue")

	Lmin = np.min([np.min(polymer[:,0]-1), np.min(polymer[:,1]-1), np.min(polymer[:,2]-1)])
	Lmax = np.min([np.max(polymer[:,0]+1), np.max(polymer[:,1]+1), np.max(polymer[:,2]+1)])
	
	ax.set_xlim(Lmin, Lmax)
	ax.set_ylim(Lmin, Lmax)
	ax.set_zlim(Lmin, Lmax)

	# save the image
	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
	"""

	# end the clock
	end = time.time()
	print(f"Time for computation is {end-start} seconds.", flush=True)
