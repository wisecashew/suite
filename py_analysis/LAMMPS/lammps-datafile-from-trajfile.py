#!/home/satyend/.conda/envs/phase/bin/python

import re
import numpy as np
import pandas as pd
import shutil
import argparse

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument ("--trajfile", dest='traj', action='store', type=str, help="Name of trajfile.")
parser.add_argument ("--datafile", dest='data', action='store', type=str, help="Name of original datafile.")
parser.add_argument ("--new-datafile", dest='ndf', action='store', type=str, help="Name of new datafile.")
parser.add_argument ("--timestep", dest='ts', action='store', type=int, help="timestep to extract from trajfile.")

args = parser.parse_args()

if __name__ == "__main__":

	data = dict()
	data["id"]     = []
	data["type"]   = []
	data["x"]      = []
	data["y"]      = []
	data["z"]      = []

	check_stepnum         = False
	before_extract_coords = False
	extract_coords        = False

	print (f"Opening the trajectory file...", flush=True)
	traj = open (args.traj, 'r')

	idx = -1
	for line in traj:
		idx+=1
		if (not extract_coords) and  re.match ("ITEM: TIMESTEP", line):
			check_stepnum = True
			continue

		elif (not check_stepnum):
			continue

		elif extract_coords and re.match("ITEM: TIMESTEP", line):
			break

		elif (not before_extract_coords) and check_stepnum:
			# print (line)
			if re.match("^"+str(args.ts)+"$", line):
				stepnum = int (line[0])
				before_extract_coords = True
				# print(f"line = {line}")
				continue
			else:
				check_stepnum = False
				continue

		elif check_stepnum and before_extract_coords and re.match(r"ITEM: ATOMS id type x y z", line):
			extract_coords = True
			continue

		elif extract_coords:
			# print (line)
			L = line.strip().split()
			data["id"].append(int(L[0]))
			data["type"].append(int(L[1]))
			data["x"].append(float(L[2]))
			data["y"].append(float(L[3]))
			data["z"].append(float(L[4]))
			continue

		else:
			continue

	traj.close ()

	if extract_coords:
		print (f"Was the timestep found? {extract_coords}.", flush=True)
	else:
		print (f"Was the timestep found? {extract_coords}.", flush=True)
		exit ()



	datafile = open(args.data, 'r')
	new_df   = open(args.ndf, 'w')

	atom_loc = False

	idx = -1
	c = 0
	for line in datafile:
		# print(line)
		idx+=1
		if re.match("Atoms", line):
			atom_loc = True
			new_df.write(line)
			continue

		elif re.match("Bonds", line):
			atom_loc = False
			new_df.write(line)
			continue

		elif line.strip() == "":
			new_df.write(line)
			continue

		elif (not atom_loc):
			new_df.write(line)
			continue

		elif atom_loc:
			# print (line)
			mol_id = int((line.strip().split())[1])
			q      = float((line.strip().split())[3])
			new_df.write("{:<7} {:<7} {:<7} {: >1.6f} {: <3.6f} {: <3.6f} {: <3.6f} {:<1} {:<1} {:<1}\n".\
			format(data["id"][c], mol_id, data["type"][c], q, data["x"][c], data["y"][c], data["z"][c], 0, 0, 0))
			c+=1
			continue

	datafile.close()
	new_df.close()

