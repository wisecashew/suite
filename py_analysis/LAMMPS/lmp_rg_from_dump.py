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

	my_settings = lmp_settings.Settings(settings)
	my_settings.read_settings()

	my_data     = lmp_data.Data(topo)
	my_data.read_data()

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	print("Begin selecting some atoms for the polymer...", end=' ', flush=True)
	polymer = u.select_atoms('resid 1')
	print("done!")
	Rg_mda      = list()
	Rg_manual   = list()
	mass_mda    = 0
	mass_manual = 0

	for atom in polymer.atoms:
		mass_mda    += atom.mass
		mass_manual += my_data.atoms[int(atom.type)]["mass"]

	print(f"mass_manual = {mass_manual}")
	print(f"mass_mda    = {mass_mda}")

	for idx, ts in enumerate(u.trajectory):
		rcom = np.array([0.0, 0.0, 0.0])
		polymer.unwrap()
		# positions = polymer.positions
		# for atom in polymer.atoms:
		# 	rcom += atom.mass * atom.position
		# rcom = rcom/mass_manual
		# Rg = 0
		# for atom in polymer.atoms:
		# 	Rg += atom.mass * np.linalg.norm(atom.position - rcom) ** 2
		# Rg = Rg/mass_manual
		# Rg = np.sqrt(Rg)
		# Rg_manual.append(Rg)
		Rg_mda.append(polymer.radius_of_gyration())
		print(f"ts: {ts}, {polymer.radius_of_gyration()}")

	print(f"mda    = {Rg_mda}")
	print(f"manual = {Rg_manual}")
