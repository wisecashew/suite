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
parser.add_argument("--latt-pkl", dest="lpkl", type=str, action="store", help="enter address of processed lattice.", required=True)
parser.add_argument("--dump",     dest="dump", type=str, action="store", help="dump of contacts.",                   required=True)
args = parser.parse_args()

if __name__=="__main__":

	# start the clock
	start = time.time()

	# extract the lattice object
	with open(args.lpkl, "rb") as latt_file:
		latt = pickle.load(latt_file)

	# set up the dump file
	dumpfile = open(args.dump, 'w')
	dumpfile.write("timestep Nmm Nms Nmsa Ns\n")

	# calculate contacts
	for idx,ts in enumerate(latt):
		print(f"Timestep = {ts}...", flush=True)
		timestep = int(ts*1e+4)
		Nmm  = latt[ts]["Nmm"]
		Nms  = 26*40 - 2*Nmm
		Nmsa = latt[ts]["Nmsa"]
		all_neighbors    = latt[ts]["polymer"][:,np.newaxis,:] + ll.shifts[np.newaxis,:,:]
		all_neighbors    = all_neighbors.reshape(-1,3)
		unique_neighbors = np.unique(all_neighbors, axis=0)
		mask             = ~np.any(np.all(unique_neighbors[:,np.newaxis,:] == latt[ts]["polymer"], axis=2), axis=1)
		unique_neighbors = unique_neighbors[mask]
		Ns               = len(unique_neighbors)

		if idx == len(latt)-1:
			dumpfile.write(f"{int(ts*1e+4)} {Nmm} {Nms} {Nmsa} {Ns}")
		else:
			dumpfile.write(f"{int(ts*1e+4)} {Nmm} {Nms} {Nmsa} {Ns}\n")
	dumpfile.close()

	# end the clock
	end = time.time()
	print(f"Time for computation is {end-start} seconds.", flush=True)
