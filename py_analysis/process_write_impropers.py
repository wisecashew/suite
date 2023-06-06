#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np
import re
import time 

# set up a dictionary that goes from id to atom 
def set_up_atom_map (lmps2type):
	atoms = dict()
	f = open (lmps2type, 'r')
	for line in f:
		pattern_first_column  = r'^(\S+)'
		pattern_second_column = r'\S+\s+(\S+)'
		key = re.findall (pattern_first_column, line)
		# print (f"key = {key}")
		value = re.findall (pattern_second_column, line)
		atoms[key[0]] = value[0]

	# print (atoms)
	return atoms

# now go to the data file and get all the impropers of type C82, C177, N180, O178
# create a data and settings file parser 

def data_settings_parser (data_file, settings_file, lmps2type):

	atom_map  = set_up_atom_map (lmps2type)
	# print (atom_map)
	molecules = dict()
	bonds     = set()

	f = open (data_file, 'r')
	parse_flag  = False
	atom_flag   = False
	bond_flag   = False
	num_pattern = r'^[0-9]+'
	pattern_first_column  = r'^\s*(\d+)\s+'
	pattern_second_column = r'^\s*(?:\S+\s+)(\d+)\s+'
	pattern_third_column  = r'^\s*(?:\S+\s+){2}(\d+)\s+'
	pattern_fourth_column = r'^(?:\S+\s+){3}(\d+)'

	for line in f:
		if "Atoms" in line:
			parse_flag = True
			atom_flag  = True

		elif "Bonds" in line:
			atom_flag = False
			bond_flag = True

		elif "Angles" in line:
			break

		if parse_flag:
			if atom_flag:
				if re.match (num_pattern, line):
					# get molecule number 
					key_molnum = re.findall (pattern_second_column, line)
					if molecules.get(key_molnum[0]) is None:
						molecules[key_molnum[0]] = dict()

					# get serial numbers 
					key_srnum = re.findall (pattern_first_column, line)
					molecules[key_molnum[0]][key_srnum[0]] = None

					# get atom type 
					key_atomtype = re.findall (pattern_third_column, line)
					# print (f"line = {line}")
					# print (key_atomtype)
					molecules[key_molnum[0]][key_srnum[0]] = (key_atomtype[0], atom_map[key_atomtype[0]])

			if bond_flag:
				if re.match (num_pattern, line):
					# get two atoms of the polymers 
					key_p1 = re.findall (pattern_third_column , line)
					key_p2 = re.findall (pattern_fourth_column, line)
					key_tup = tuple([key_p1[0], key_p2[0]])
					bonds.add (key_tup)

		else:
			continue

	return atom_map, molecules, bonds


def improper_finder (data_file, settings_file, lmps2type, num_particles_per_monomer):

	atom_map, molecules, bonds = data_settings_parser (data_file, settings_file, lmps2type)

	# find all the polymers
	polymers = dict()
	for key in molecules:
		if len (molecules[key].keys()) > 100:
			polymers[key] = dict()
			polymers[key]["N"] = round((len(molecules[key].keys())-2)/num_particles_per_monomer)
			polymers[key]["particles"] = molecules[key]

	# find the improper dihedrals
	impr_dihedrals = set()
	initial = 0
	for key in polymers:
		print (f"polymer number = {key} with degree of polymerization = " +str(polymers[key]["N"]))
		for monomer_num in range (polymers[key]["N"]):
			if monomer_num == 0 or monomer_num == polymers[key]["N"]-1:
				tot_part = num_particles_per_monomer + 1
				idih = []
				for j in range (tot_part):
					if monomer_num == polymers[key]["N"]-1:
						if polymers[key]["particles"][str(initial+j+1)][1] in ["O178", "N180", "C177", "C81"]:
							idih.append (str(initial+j+1))
					else:
						if polymers[key]["particles"][str(initial+j+1)][1] in ["O178", "N180", "C177", "C82"]:
							idih.append (str(initial+j+1))

				if monomer_num == polymers[key]["N"]-1:
					# print (idih)
					# find the sr num corresponding to C177
					c81_list = []
					for sr_num in idih:
						if polymers[key]["particles"][sr_num][1] == "C177":
							c177 = sr_num
						elif polymers[key]["particles"][sr_num][1] == "C81":
							c81_list.append (sr_num)

					# check which C81 is bonded to C177
					for c81 in c81_list:
						check = tuple (sorted([c81, c177]))
						if check not in bonds:
							idih.remove(c81)

				if len (idih) != 4:
					print ("Something wrong with improper dihedrals for monomer_num = " + str(monomer_num))
					print (idih)
					exit()

				# sort the improper from C82/C81 to C177 to N180 to O178
				sorted_idih = [None, None, None, None]
				for i in range(4):
					if polymers[key]["particles"][idih[i]][1] == "N180":
						sorted_idih[2] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "O178":
						sorted_idih[3] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C177":
						sorted_idih[1] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C82":
						sorted_idih[0] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C81" and monomer_num == polymers[key]["N"]-1:
						sorted_idih[0] = idih[i]


				impr_dihedrals.add (tuple(sorted_idih))
				initial += tot_part 

			else:
				tot_part = num_particles_per_monomer
				idih = []
				for j in range (tot_part):
					if polymers[key]["particles"][str(initial+j+1)][1] in ["O178", "N180", "C177", "C82"]:
						idih.append (str(initial+j+1))

				if len (idih) != 4:
					print ("Something wrong with improper dihedrals for monomer_num = " + str(monomer_num))
					print (idih)
					exit()

				# sort the improper from C82/C81 to C177 to N180 to O178
				sorted_idih = [None, None, None, None]
				for i in range(4):
					if polymers[key]["particles"][idih[i]][1] == "N180":
						sorted_idih[2] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "O178":
						sorted_idih[3] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C177":
						sorted_idih[1] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C82":
						sorted_idih[0] = idih[i]
					elif polymers[key]["particles"][idih[i]][1] == "C81" and monomer_num == polymers[key]["N"]-1:
						sorted_idih[0] = idih[i]

				impr_dihedrals.add (tuple(sorted_idih))
				initial += tot_part

			# print ("moving to next monomer unit = " + str(monomer_num) + ", initial = " + str(initial))

	return impr_dihedrals


def header_index_finder_datafile (filename):

	f = open (filename, 'r')
	for idx, line in enumerate(f):
		if re.findall (r'dihedral types', line):
			return idx

	return None

def header_inserter_datafile (filename, line_index, improper_num):

	f = open (filename, 'r')
	content = f.readlines ()
	content.insert (line_index, f"{improper_num} impropers\n1 improper types\n")
	f.close ()
	with open (filename+".new", 'w') as f:
		content = "".join(content)
		f.write (content)
	return

def tail_improper_inserter_datafile (filename, impropers):

	f = open (filename+".new", 'a')
	f.write ("\n")
	f.write ("Impropers\n")
	f.write ("\n")
	for idx, element in enumerate(impropers):
		f.write (f"{idx+1}    1    {element[0]}    {element[1]}    {element[2]}    {element[3]}\n")
	f.write ("\n")
	f.close ()
	return

def tail_improper_inserter_settingsfile (filename):
	f = open (filename, 'a')
	f.write  ("\n")
	f.write  ("# impropers \n")
	f.write  ("improper_coeff  1   21  0\n")
	f.write  ("improper_coeff  1   aa  0   0   0   120 120 120")
	f.close  ()
	return 

import argparse
parser = argparse.ArgumentParser(description="Reads the data file and settings file created by gen_data for PNIPAM and then adds impropers to those files, externally.")
parser.add_argument ('--data', dest='d', action='store', type=str, help="Name of datafile to edit.")
parser.add_argument ('--settings', dest='s', action='store', type=str, help="Name of settings file to edit.")
parser.add_argument ('--map', dest='m', action='store', type=str, help="Name of mapping file.")
parser.add_argument ('-n', dest='n', action='store', type=int, help="Number of particles per monomer unit.")

args = parser.parse_args ()


if __name__=="__main__":

	start = time.time()
	data_file = args.d      # "sys.data"
	settings_file = args.s  # "sys.settings"
	lmps2type = args.m      # "lmps2type.map"
	num_particles_per_monomer = args.n
	impr_dihedral_set = improper_finder (data_file, settings_file, lmps2type, num_particles_per_monomer)

	# find the line index where the improper types need to be plugged in 
	# this happens write after "dihedral types"

	line_index = header_index_finder_datafile (args.d)
	header_inserter_datafile (args.d, line_index+1, len(impr_dihedral_set))
	tail_improper_inserter_datafile (args.d, impr_dihedral_set)
	tail_improper_inserter_settingsfile (args.s)

	stop = time.time()
	print (f"Time required for computation is {stop-start} seconds.")
