#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
import lmp_data
import lmp_settings
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Print out the radius of gyration.")
parser.add_argument("--datafile",      dest="df",    type=str, action="store", required=True, help="enter address of datafile.")
parser.add_argument("--settingsfile",  dest="sf",    type=str, action="store", required=True, help="enter address of settingsfile.")
parser.add_argument("--trajfile",      dest="tf",    type=str, action="store", required=True, help="enter address of coordsfile.", metavar="coords.lammpstrj")
parser.add_argument("--rg-dump",       dest="rg",    type=str, action="store", required=True, help="enter dump file for rg.")
args = parser.parse_args()

if __name__=="__main__":
	
	settings = args.sf
	topo     = args.df
	traj     = args.tf
	dt       = 2

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	print("Begin selecting some atoms for the polymer...", end=' ', flush=True)
	polymer = u.select_atoms('resid 1')

	rg_dumpfile = open(args.rg, 'w')
	for idx, ts in enumerate(u.trajectory):
		polymer.unwrap()
		Rg_mda = polymer.radius_of_gyration()
		rg_dumpfile.write(f"{ts.time/2} {Rg_mda}\n")
	rg_dumpfile.close()

