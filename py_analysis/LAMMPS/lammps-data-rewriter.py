#!/home/satyend/.conda/envs/phase/bin/python

import re
import numpy as np
import pandas as pd
import shutil
import argparse

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument ("--datafile", dest='data', action='store', type=str, help="Name of datafile.")
parser.add_argument ("--edited", dest='edit', action='store', type=str, help="Name of new, edited file.")
parser.add_argument ("--x", dest='x', action='store', type=float, help="Length along X.")
parser.add_argument ("--y", dest='y', action='store', type=float, help="Length along Y.")
parser.add_argument ("--z", dest='z', action='store', type=float, help="Length along Z.")

args = parser.parse_args()


def get_atom_bonds_map(datafile):

	data = dict()
	data["srno"]   = []
	data["molid"]  = []
	data["atomid"] = []
	data["q"]      = []
	data["x"]      = []
	data["y"]      = []
	data["z"]      = []
	data["bx"]     = []
	data["by"]     = []
	data["bz"]     = []

	bonds = dict()

	atoms_index = False
	bonds_index = False
	count = 0
	f = open (args.data, 'r')
	for idx,line in enumerate(f):
		if re.match ("Atoms\n", line):
			atoms_index = True
			continue
		elif re.match ("Bonds\n", line):
			# print (line)
			bonds_index = True
			atoms_index = False
			continue

		elif re.match ("Angles", line):
			break

		if re.match (r'\s*\S+(\s+\S+){9}\s*$', line) and atoms_index:
			contents = line.strip().split()
			contents = [ float(elem) for elem in contents ]
			data["srno"].append  (int(contents[0]))
			data["molid"].append (int(contents[1]))
			data["atomid"].append(int(contents[2]))
			data["q"].append(contents[3])
			data["x"].append(contents[4])
			data["y"].append(contents[5])
			data["z"].append(contents[6])
			data["bx"].append(contents[7])
			data["by"].append(contents[8])
			data["bz"].append(contents[9])
			count += 1
			continue

		elif re.match(r'\d+\s+\d+\s+\d+\s+\d+', line) and bonds_index:
			contents = line.strip().split()
			atom_nums = [int(contents[2]), int(contents[3])]
			if atom_nums[0] in bonds:
				bonds[atom_nums[0]].append(atom_nums[1])
			else:
				bonds[atom_nums[0]] = [atom_nums[1]]

			if atom_nums[1] in bonds:
				bonds[atom_nums[1]].append(atom_nums[0])
			else:
				bonds[atom_nums[1]] = [atom_nums[0]]

	f.close()

	return data, bonds

def unwrap(polymer_wrapped, polymer_unwrapped, data, bond_tree, idx, x, y, z):

	# start with the first index
	if idx == 0:
		polymer_unwrapped[idx] = polymer_wrapped[idx]

	bonded = bond_tree[data["srno"][idx]]
	print(idx)
	print(bonded)
	for bonded_atom in bonded:
		print(bond_tree[data["srno"][bonded_atom]])
		bond_tree[data["srno"][bonded_atom]].remove(idx)
		displacement = polymer_wrapped[bonded_atom-1] - polymer_wrapped[idx]
		disp = np.array([0,0,0])

		# Calculate displacement vector using vectorized operations
		mask = abs(displacement) > np.array([x/2, y/2, z/2])
		disp[mask & (displacement > 0)] = displacement[mask & (displacement > 0)] - np.array([x, y, z])
		disp[mask & (displacement < 0)] = displacement[mask & (displacement < 0)] + np.array([x, y, z])

		polymer_unwrapped[bonded_atom - 1] = polymer_unwrapped[idx] + disp

	del bond_tree[idx]


	if len(bond_tree) == 0:
		return polymer_unwrapped
	else:
		for bonded_atom in bonded:
			unwrap(polymer_wrapped, polymer_unwrapped, data, bond_tree, bonded_atom-1, x, y, z)

	return polymer_unwrapped 


if __name__ == "__main__":

	datafile = args.data

	data, bond_tree = get_atom_bonds_map(datafile)

	# print(bond_tree)

	x = 75
	y = 75
	z = 75

	# now that we have the bonds, we can now start unwrapping molecules
	# 1. count the number of molecules in the system
	# only choose the first molecule
	polymer_atoms_index = [i for i, x in enumerate(data["molid"]) if x == 1]
	polymer_atoms = [data["srno"][i] for i in polymer_atoms_index]
	
	polymer_wrapped = np.zeros((len(polymer_atoms), 3))
	for p_srno in polymer_atoms:
		polymer_wrapped[p_srno-1] = np.array([data["x"][p_srno-1], data["y"][p_srno-1], data["z"][p_srno-1]])

	# 2. unwrap the polymer
	start_idx = 0
	polymer_unwrapped = np.zeros(polymer_wrapped.shape)
	polymer_unwrapped = unwrap(polymer_wrapped, polymer_unwrapped, data, bond_tree, start_idx, x, y, z)

	for idx,loc in enumerate(polymer_unwrapped):
		data["x"][idx] = polymer_unwrapped[idx][0]
		data["y"][idx] = polymer_unwrapped[idx][1]
		data["z"][idx] = polymer_unwrapped[idx][2]


	df = pd.DataFrame.from_dict (data)

	# get center of mass of polymer
	df = df[df["molid"]==1]
	x_com = np.mean (df["x"].values)
	y_com = np.mean (df["y"].values)
	z_com = np.mean (df["z"].values)
	r_com = np.array ([x_com, y_com, z_com])
	center = np.array([args.x/2, args.y/2, args.z/2])

	translate = center - r_com 

	df["x"] = df["x"] + translate[0]
	df["y"] = df["y"] + translate[1]
	df["z"] = df["z"] + translate[2]

	count = 0
	f = open (args.data, 'r')
	g = open (args.edit, 'w')
	atom_index = False
	for line in f:
		
		if re.match ("Atoms\n", line):
			atom_index = True
		elif re.match ("Bonds\n", line):
			atom_index = False
		
		if re.match (r'^\s*\S+(\s+\S+){9}\s*$', line) and atom_index and count < 572:
			g.write ("{} {} {} {} {} {} {} {} {} {}\n".format(df["srno"].values[count], df["molid"].values[count], df["atomid"].values[count], df["q"].values[count], df["x"].values[count], df["y"].values[count], df["z"].values[count], df["bx"].values[count], df["by"].values[count], df["bz"].values[count]))
			count += 1
		else:
			g.write (line)

	f.close()
	g.close()

	shutil.copyfile (args.edit, args.data)

