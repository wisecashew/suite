#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
import re
import argparse

parser = argparse.ArgumentParser(description="Get the information necessary for persistence length.")
parser.add_argument("--settingsfile",  dest="sf",    type=str, action="store", required=True, help="enter address of settingsfile.")
parser.add_argument("--datafile",      dest="df",    type=str, action="store", required=True, help="enter address of datafile.")
args = parser.parse_args()

if __name__=="__main__":
	
	# --- set up the inputs --- 
	settings = args.sf       # this is the settings file
	dt       = 2             # this is the size of the timestep
	topo     = args.df       # this is the data file

	# --- get the trajectory ---
	traj     = f"POLYMER_COORD_FILES/polymer.0.lammpstrj"	# get the trajectory file

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	# ~~~ this is my particular case ~~~
	backbone_selection = "id " # this is the start of the backbone
	for i in range(40):
		if i == 0:
			backbone_selection += f"{19*i + 1} {19*i + 5} "
		elif i == 39:
			backbone_selection += f"{19*i + 2} {19*i + 5}"
		else:
			backbone_selection += f"{19*i + 2} {19*i + 5} "
	print(f"\nbackbone_selection = {backbone_selection}") # printing out the backbone
	backbone = u.select_atoms(backbone_selection)

	for atom in backbone[1::2]:
		print(f"sr. no. = {atom.index+1}, type = {atom.type}.")
