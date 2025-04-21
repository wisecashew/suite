#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import nglview as nv
import argparse
import pickle
import MDAnalysis.transformations as xforms
from mpl_toolkits.mplot3d import Axes3D
import latt_lmp as ll

parser = argparse.ArgumentParser(description="Compute hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--datafile",     dest="df",   type=str, action="store", help="enter address of datafile.",                     required=True)
parser.add_argument("--settingsfile", dest="sf",   type=str, action="store", help="enter address of settingsfile.",                 required=True)
parser.add_argument("--coords",       dest='c',    type=str, action="store", help="enter address of coordinate file.",              required=True)
parser.add_argument("--Lpkl",         dest="lpkl", type=str, action="store", help="enter address of final processed lattice.",      required=True)
args = parser.parse_args()

if __name__=="__main__":

	# start the clock
	start = time.time()

	# get the names of the files
	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	# generate the universe object
	print(f"Making the universe...", end=' ', flush=True)
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"made!", flush=True)

	# unwrap the trajectory
	polymer = u.select_atoms("resid 1")
	other   = u.select_atoms("not resid 1")
	_all    = u.select_atoms("all")

	# set up the transforms
	transforms = [mda.transformations.unwrap(_all), mda.transformations.center_in_box(polymer, wrap=False), mda.transformations.wrap(_all)]
	u.trajectory.add_transformations(*transforms)

	# define the lattice_lammps object
	LL = ll.Lattice_LMP(u, dict())
	LL.latticify_dry()

	# calculate contacts
	for ts in LL.lattice:
		print(f"Calculating contacts @ timestep = {ts}. Calculating Nmm...", end=' ', flush=True)
		LL.get_Nmm (ts)
		print(f"done!", end=' ', flush=True)


	my_lattice = open(args.lpkl, "wb")
	pickle.dump(LL, my_lattice)
	my_lattice.close()

	end = time.time()
	print(f"Time for computation is {end-start} seconds.", flush=True)
