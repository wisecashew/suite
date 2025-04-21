#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
import MDAnalysis.transformations as xforms
import re
import lmp_data
import lmp_settings
import numpy as np
import pandas as pd
import argparse
import pickle
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


parser = argparse.ArgumentParser(description="Get the information necessary for persistence length.")
parser.add_argument("--settingsfile",  dest="sf",    type=str, action="store", required=True, help="enter address of settingsfile.")
parser.add_argument("--datafile",      dest="df",    type=str, action="store", required=True, help="enter address of datafile.")
parser.add_argument("--pl-dump",       dest="pl",    type=str, action="store", required=True, help="enter dump file for persistence length.")
parser.add_argument("--sim",           dest="sim",   type=int, action="store", required=True, help="enter index of simulation.")
args = parser.parse_args()

def get_Tidx(df, time_point, replica_id):
	# Find the row index where time matches
	row_idx = np.where(df.iloc[:, 0].values == time_point)[0][0]

	# Find the column index of the temperature index
	col_idx = np.where(df.iloc[row_idx, 1:].values == replica_id)[0][0]

	return col_idx


if __name__=="__main__":
	
	name_str = ["time"]
	for rep_idx in range(40):
		name_str.append(str(rep_idx))

	# get the output files
	df = pd.read_csv("TOT_OUTPUT_FILES/output.dump", engine="python", names=name_str, sep='\s+', skiprows=0)

	# --- set up the inputs --- 
	settings = args.sf       # this is the settings file
	dt       = 2             # this is the size of the timestep
	topo     = args.df       # this is the data file
	sim_idx  = args.sim      # this is the index of the replica

	# get the selection string
	print("Begin selecting some atoms for the polymer...", end=' ', flush=True)
	backbone_selection = "id " # this is the start of the backbone
	for i in range(40):
		if i == 0:
			backbone_selection += f"{19*i + 1} {19*i + 5} "
		elif i == 39:
			backbone_selection += f"{19*i + 2} {19*i + 5}"
		else:
			backbone_selection += f"{19*i + 2} {19*i + 5} "
	print(f"\nbackbone_selection = {backbone_selection}") # printing out the backbone

	# set up the datastructures
	max_delta = 39			# this is the maximum separation between monomers on a 40-mer
	pl_database = dict()	# this is a database that stores all the correlations over temperature data in a database
	T_range = np.linspace(275, 360, 40)
	for Tidx in range(40):
		pl_database[Tidx] = dict()
		pl_database[Tidx]["bond_length"] = []
		for delta in range(1,max_delta):
			pl_database[Tidx][delta] = []

	# start setting up the thing
	start = time.time() 

	# get the trajectory
	traj     = f"POLYMER_COORD_FILES/polymer.{sim_idx}.lammpstrj"	# get the trajectory file

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	backbone = u.select_atoms(backbone_selection) # get the backbone
	_all     = u.select_atoms("all")
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	# transform the backbone: this keeps the backbone in the center of the box
	transforms = [mda.transformations.unwrap(_all), mda.transformations.center_in_box(backbone, wrap=False), mda.transformations.wrap(_all)]
	u.trajectory.add_transformations(*transforms)
	
	for idx, ts in enumerate(u.trajectory):
		print(f"At timestep {ts.time}.", flush=True)
		Tidx = get_Tidx(df, int(ts.time/2), sim_idx)								# get the temperature at which this particular configuration was sampled
		backbone_elements = np.array(backbone.positions)							# get the coordinates as a numpy array
		backbone_elements = backbone_elements[1::2]									# only consider the carbons that are attached to the oxygen
		bond_elements     = backbone_elements[1:] - backbone_elements[:-1]			# getting the distances between the adjacent atoms
		norms             = np.linalg.norm(bond_elements, axis=1, keepdims=True)	# get the lengths of the bonds between atoms
		pl_database[Tidx]["bond_length"].append(np.mean(norms))						# get the average bond length for this point
		unit_vectors      = bond_elements / norms									# normalize the bond vectors

		# there are 80 atoms and 79 bonds.
		# therefore, when getting correlations, delta can go from 1 to 78.
		# at delta = 0, correlation is exactly 1.
		for delta in range(1, max_delta):
			dot_products = np.sum(unit_vectors[:-delta] * unit_vectors[delta:], axis=1)  # Vectorized dot products
			pl_database[Tidx][delta].extend(dot_products.tolist())                       # Directly update the database

	stop = time.time()
	print(f"Time for #{sim_idx} is {stop-start} seconds.", flush=True)
	
	pl_file = open(args.pl, 'wb')
	pickle.dump(pl_database, pl_file)
	pl_file.close()

