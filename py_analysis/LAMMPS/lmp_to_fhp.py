#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import latt_lmp as ll
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import numpy as np
import time
import argparse
import pickle
import MDAnalysis.transformations as trans
from mpl_toolkits.mplot3d import Axes3D
import sys
np.set_printoptions(threshold=sys.maxsize)

parser = argparse.ArgumentParser(description="Compute hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--datafile",     dest="df",    type=str, action="store", help="enter address of datafile.",        required=True)
parser.add_argument("--settingsfile", dest="sf",    type=str, action="store", help="enter address of settingsfile.",    required=True)
parser.add_argument("--coords",       dest='c',     type=str, action="store", help="enter address of coordinate file.", required=True)
parser.add_argument("--Hpkl",         dest="Hpkl",  type=str, action="store", help="enter address of hydrogen bonding dump.",     required=True)
args = parser.parse_args()

if __name__=="__main__":

	# start the clock
	start = time.time()

	# get the names of the files
	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	# define some figures
	# this is for the MD coordinates
	fig = plt.figure(num=0)
	ax  = fig.add_subplot(111, projection="3d")

	# this is for the lattice coordinates
	flatt  = plt.figure(num=1)
	axlatt = flatt.add_subplot(111, projection="3d")

	# generate the universe object
	print(f"Making the universe...", end=' ', flush=True)
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"made!", flush=True)

	# unwrap the trajectory
	polymer = u.select_atoms("resid 1")
	other   = u.select_atoms("not resid 1")
	_all    = u.select_atoms("all")

	# transform the universe
	transforms = [mda.transformations.unwrap(_all), mda.transformations.center_in_box(polymer, wrap=False), mda.transformations.wrap(_all)]
	u.trajectory.add_transformations(*transforms)

	# get the hbonding info
	file = open(args.Hpkl, "rb")
	data = pickle.load(file)
	file.close()

	hbonds = dict()
	hbonds["Op_with_water"] = np.array(data[0].results.hbonds[:,0:4], dtype=int)
	hbonds["Np_with_water"] = np.array(data[1].results.hbonds[:,0:4], dtype=int)
	hbonds["Hp_with_water"] = np.array(data[2].results.hbonds[:,0:4], dtype=int)
	hbonds["Op_with_Hp"]    = np.array(data[3].results.hbonds[:,0:4], dtype=int)

	# get the colormaps and colors
	cmap   = cm.get_cmap("tab20", 40)
	colors = [cmap(i) for i in range(40)]

	# define the containers for lattice positions and center of masses of monomers
	c_of_m = list()
	latt_m = list()
	w_mols = list()

	# get the start and end indices
	get_start = lambda m: 1  if m == 0 else 21 + 19 * (m-1)
	get_end   = lambda m: get_start(m)+19 if (m == 0 or m == 39) else get_start(m)+18

	# traverse through the trajectory
	for idx,ts in enumerate(u.trajectory):
		print(f"@ ts = {ts.time}.", flush=True)
		c_of_m.clear()
		latt_m.clear()
		w_mols.clear()
		for midx in range(40):
			sel_str    = f"bynum {get_start(midx)}:{get_end(midx)}"
			selection  = u.select_atoms(sel_str)
			c_of_m.append(selection.center_of_mass())
			if midx == 0:
				latt_m.append(np.array([0,0,0], dtype=int))
			else:
				ll.get_lat(c_of_m[-1], c_of_m[-2], latt_m)

		# get the hydrogen bonds at this timestep
		hb_Ow  = hbonds["Op_with_water"][hbonds["Op_with_water"][:,0] == idx]
		hb_Nw  = hbonds["Np_with_water"][hbonds["Np_with_water"][:,0] == idx]
		hb_net = np.hstack((hb_Ow, hb_Nw))

		# get the water molecules that are H-bonding
		for hb in hb_net:
			water = u.select_atoms(f"bynum {hb[1]} or bynum {hb[2]}")
			resid = water.residues.resids[0]
			w_mols.append(resid) 
		w_mols = list(set(w_mols))

		# loop through each hb
		for idx,hb in enumerate(hb_net):
			midx = (hb[-1]-2)//19 if (hb[-1]-2)//19 < 40 else 39
			# get the molecule number from
			water = u.select_atoms(f"resid {w_mols[idx]}")
			resid = water.residues.resids[0]
			water = u.select_atoms(f"resid {resid}")
			w_com = water.center_of_mass()
			ll.get_lat(w_com, latt_m[midx], latt_m)

		# urows, indices, counts = np.unique(np.array(latt_m), axis=0, return_index=True, return_counts=True)
		# if np.any(counts > 1):
		# 	print("Yes, duplicates.", flush=True)
		# else:
		# 	print(f"No duplicates!", flush=True)




