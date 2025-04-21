import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import argparse
import time

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the hydrogen bonds over time.")
parser.add_argument('--coords',    dest='c',     action='store', type=str,   help='Name of the coordinate dump.', default=None)
parser.add_argument('--data',      dest='data',  action='store', type=str,   help='Name of the data file.', default=None)
parser.add_argument('--image',     dest='img',   action='store', type=str,   help='Name of image to be created.', default=None)
args = parser.parse_args()

if __name__=="__main__":

	start = time.time()

	try:
		u = mda.Universe(args.data, args.c, format="LAMMPS", lengthunit="angstrom", timeunit="ps")
	except:
		u = mda.Universe(args.data, args.c, atom_style="id type x y z", lengthunit="angstrom", timeunit="ns") #, lengthunit="angstrom", timeunit="ns")

	print("what is the universe object?", flush=True)
	print(f"u = {u}")
	print()

	# traverse through the trajectory with the following:
	print("what is the trajectory object inside the universe object?")
	print(f"u.trajectory = {u.trajectory}", flush=True)
	print()

	print(f"what is this DCDReader object?")
	print(f"I don't know for sure. It is some kind of an interator object. So let's iterate through it.")
	print()

	print(f"u.trajectory.time = {u.trajectory.time}")

	for idx,ts in enumerate(u.trajectory):
		print(f"Element ts (of u.trajectory) = {ts}.")
		print(f"Printing out some properties of each timestep...")
		print(f"ts.n_atoms              = {ts.n_atoms}")
		print(f"ts.time                 = {ts.time}")
		print(f"ts.dt                   = {ts.dt}")
		print(f"ts.has_positions        = {ts.has_positions}")
		print(f"ts.has_velocities       = {ts.has_velocities}")
		print(f"ts.has_forces           = {ts.has_forces}")
		print(f"ts.dimensions           = {ts.dimensions}")
		print()
		print("========================")
		print()

	print(f"Check what u.dimensions outputs.")
	print(f"u.dimensions = {u.dimensions}")
	print()
	print("========================")
	print()

	print("Check what the atoms class is in u.")
	print(f"u.atoms = {u.atoms}")

	'''
	for idx, at in enumerate(u.atoms):
		print(f"at = {at}")
		print(f"\tat.bonded_atoms={at.bonded_atoms}, \n\
			at.type = {at.type}, at.resid = {at.resid}, at.index = {at.index}, at.mass = {at.mass}, at.charge = {at.charge}, \n\
			at.position = {at.position},\n \
			at.get_connections = {at.get_connections}.")
	'''
	# time to play with molecules in the universe
	print("Begin selecting some atoms for the polymer...")
	polymer = u.select_atoms('resid 1')
	print(f"polymer = {polymer}")
	Rg = []

	for ts in u.trajectory:
		polymer.unwrap()
		Rg.append(polymer.radius_of_gyration(wrap=False))

	print(f"Rg = {Rg}")

	# test = u.select_atoms('all')
	# test = test[test.resids > 1]

	# find all water molecules around the polymer
	surrounding_set = set()
	nsurr_water = []
	for ts in u.trajectory:
		surrounding_water = u.select_atoms("around 4 resid 1")
		surrounding_water = surrounding_water[surrounding_water.resids > 1]
		for atoms in surrounding_water:
			surrounding_set.add(atoms.resid)
		nsurr_water.append(len(surrounding_set))
		surrounding_set.clear()

	print(nsurr_water)

	stop = time.time()
	print(f"Time for exploration is {stop-start} seconds.", flush=True)


