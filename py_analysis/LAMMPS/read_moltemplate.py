#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import argparse
import pickle 
parser = argparse.ArgumentParser(description="Read moltemplate output.")
parser.add_argument("--datafile",     dest='df', type=str, action='store', help="enter address of datafile.")
parser.add_argument("--settingsfile", dest='sf', type=str, action='store', help="enter address of settingsfile.")
parser.add_argument("--dpkl", dest='dpkl', type=str, action='store', help="enter address of data object dump.", default=None)

args = parser.parse_args()

class Data:
	def __init__(self):
		self.natoms     = 0
		self.nbonds     = 0
		self.nangles    = 0
		self.ndihedrals = 0
		self.nimpropers = 0
		self.natom_types  = 0
		self.nbond_types  = 0
		self.nangle_types = 0
		self.ndih_types   = 0
		self.nimpr_types  = 0
		self.xdim   = [0, 0]
		self.ydim   = [0, 0]
		self.zdim   = [0, 0]
		self.atoms  = dict()
		self.bonds  = dict()
		self.angles = dict()
		self.dih    = dict()
		self.impr   = dict()
		self.atoms_serial  = dict()
		self.bonds_serial  = dict()
		self.angles_serial = dict()
		self.dih_serial    = dict()
		self.impr_serial   = dict()
		self.b_description = False
		self.b_masses      = False
		self.b_atoms       = False
		self.b_bonds       = False
		self.b_angles      = False
		self.b_dihedrals   = False
		self.b_impropers   = False
		return

	def check_section(self, line):
		if re.findall(r"^LAMMPS Description", line):
			self.b_description = True
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = False
		elif re.findall(r"^Masses", line):
			self.b_description = False
			self.b_masses      = True
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = False
		elif re.findall(r"^Atoms", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = True
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = False
		elif re.findall(r"^Bonds", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = True
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = False
		elif re.findall(r"^Angles", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = True
			self.b_dihedrals   = False
			self.b_impropers   = False
		elif re.findall(r"^Dihedrals", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = True
			self.b_impropers   = False
		elif re.findall(r"^Impropers", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = True
		return

		return
	
	def read_description(self, line):
		if re.findall("atoms", line):
			line = line.strip().split()
			self.natoms = int(line[0])
		elif re.findall("bonds", line):
			line = line.strip().split()
			self.nbonds = int(line[0])
		elif re.findall("angles", line):
			line = line.strip().split()
			self.nangles = int(line[0])
		elif re.findall("dihedrals", line):
			line = line.strip().split()
			self.ndihedrals = int(line[0])
		elif re.findall("impropers", line):
			line = line.strip().split()
			self.nimpropers = int(line[0])
		elif re.findall("atom types", line):
			line = line.strip().split()
			self.natom_types = int(line[0])
		elif re.findall("bond types", line):
			line = line.strip().split()
			self.nbond_types = int(line[0])
		elif re.findall("angle types", line):
			line = line.strip().split()
			self.nangle_types = int(line[0])
		elif re.findall("dihedral types", line):
			line = line.strip().split()
			self.ndih_types = int(line[0])
		elif re.findall("improper types", line):
			line = line.strip().split()
			self.nimpr_types = int(line[0])
		elif re.findall("xlo xhi", line):
			line = line.strip().split()
			self.xdim[0] = line[0]
			self.xdim[1] = line[1]
		elif re.findall("ylo yhi", line):
			line = line.strip().split()
			self.ydim[0] = line[0]
			self.ydim[1] = line[1]
		elif re.findall("zlo zhi", line):
			line = line.strip().split()
			self.zdim[0] = line[0]
			self.zdim[1] = line[1]
		return

	def read_masses(self, line):
		line = line.strip().split()
		self.atoms[int(line[0])] = dict()
		self.atoms[int(line[0])]["type"]   = line[-1]
		self.atoms[int(line[0])]["mass"]   = line[1]
		self.atoms[int(line[0])]["sr_nos"] = list()
		return

	def read_atoms(self, line):
		line = line.strip().split()
		self.atoms[int(line[2])]["sr_nos"].append(int(line[0]))
		self.atoms_serial[int(line[0])] = int(line[2])
		return

	def read_bonds(self, line):
		line    = line.strip().split()
		atype_1 = self.atoms_serial[int(line[2])]
		atype_2 = self.atoms_serial[int(line[3])]
		btype   = int(line[1])
		if btype in self.bonds:
			pass
		else:
			self.bonds[btype] = dict()
			self.bonds[btype]["sr_nos"]  = list()
			self.bonds[btype]["type"]    = set()
		self.bonds[btype]["type"].add((self.atoms[atype_1]["type"], self.atoms[atype_2]["type"]))
		self.bonds[btype]["sr_nos"].append(int(line[0]))
		self.bonds_serial[int(line[0])]  = btype
		return

	def read_angles(self, line):
		line   = line.strip().split()
		atype_1 = self.atoms_serial[int(line[2])]
		atype_2 = self.atoms_serial[int(line[3])]
		atype_3 = self.atoms_serial[int(line[4])]
		angtype = int(line[1])
		if angtype in self.angles:
			pass
		else:
			self.angles[angtype] = dict()
			self.angles[angtype]["sr_nos"] = list()
			self.angles[angtype]["type"]   = set()
		self.angles[angtype]["type"].add((self.atoms[atype_1]["type"], self.atoms[atype_2]["type"], self.atoms[atype_3]["type"]))
		self.angles[angtype]["sr_nos"].append(int(line[0]))
		self.angles_serial[int(line[0])] = angtype
		return

	def read_dihedrals(self, line):
		line = line.strip().split()
		atype_1 = self.atoms_serial[int(line[2])]
		atype_2 = self.atoms_serial[int(line[3])]
		atype_3 = self.atoms_serial[int(line[4])]
		atype_4 = self.atoms_serial[int(line[5])]
		dtype   = int(line[1])
		if dtype in self.dih:
			pass
		else:
			self.dih[dtype] = dict()
			self.dih[dtype]["sr_nos"] = list()
			self.dih[dtype]["type"]   = set()
		self.dih[dtype]["type"].add((self.atoms[atype_1]["type"], self.atoms[atype_2]["type"], self.atoms[atype_3]["type"], self.atoms[atype_4]["type"]))
		self.dih[dtype]["sr_nos"].append(int(line[0]))
		self.dih_serial[int(line[0])] = dtype
		return

	def read_impropers(self, line):
		line = line.strip().split()
		atype_1 = self.atoms_serial[int(line[2])]
		atype_2 = self.atoms_serial[int(line[3])]
		atype_3 = self.atoms_serial[int(line[4])]
		atype_4 = self.atoms_serial[int(line[5])]
		itype  = int(line[1])
		if itype in self.impr:
			pass
		else:
			self.impr[itype] = dict()
			self.impr[itype]["sr_nos"] = list()
			self.impr[itype]["type"]   = set()
		self.impr[itype]["type"].add((self.atoms[atype_1]["type"], self.atoms[atype_2]["type"], self.atoms[atype_3]["type"], self.atoms[atype_4]["type"]))
		self.impr[itype]["sr_nos"].append(int(line[0]))
		self.impr_serial[int(line[0])] = itype
		return

if __name__=="__main__":
	
	my_data = Data()

	with open(args.df, 'r') as f:
		for line in f:
			if line.strip() == "":
				continue
			else:
				my_data.check_section(line)
				line_split = line.strip().split()
				if re.findall(r"[a-zA-Z]+", line_split[0]):
					continue
				elif my_data.b_description:
					my_data.read_description(line)
				elif my_data.b_masses:
					my_data.read_masses(line)
				elif my_data.b_atoms:
					my_data.read_atoms(line)
				elif my_data.b_bonds:
					# print(my_data.atoms)
					my_data.read_bonds(line)
				elif my_data.b_angles:
					my_data.read_angles(line)
				elif my_data.b_dihedrals:
					my_data.read_dihedrals(line)
				elif my_data.b_impropers:
					my_data.read_impropers(line)

	with open(args.sf, 'r') as f:
		for line in f:
			if line == "":
				continue
			else:
				if re.findall(r"^pair_coeff", line):
					line = line.strip().split()
					type_ = int(line[1])
					my_data.atoms[type_]["eps"] = float(line[-2])
					my_data.atoms[type_]["sig"] = float(line[-1])
				elif re.findall(r"^bond_coeff", line):
					line = line.strip().split()
					type_ = int(line[1])
					my_data.bonds[type_]["kd"] = float(line[-2])
					my_data.bonds[type_]["r0"] = float(line[-1])
				elif re.findall(r"^angle_coeff", line):
					line = line.strip().split()
					type_ = int(line[1])
					my_data.angles[type_]["k0"] = float(line[-2])
					my_data.angles[type_]["t0"] = float(line[-1])
				elif re.findall(r"^dihedral_coeff", line):
					line = line.strip().split()
					type_ = int(line[1])
					my_data.dih[type_]["V0"] = float(line[-4])
					my_data.dih[type_]["V1"] = float(line[-3])
					my_data.dih[type_]["V2"] = float(line[-2])
					my_data.dih[type_]["V3"] = float(line[-1])
				elif re.findall(r"^improper_coeff", line):
					line = line.strip().split()
					type_ = int(line[1])
					my_data.impr[type_]["k"] = float(line[-2])
					my_data.impr[type_]["phi"] = float(line[-1])
				elif re.findall(r"charge", line):
					line = line.strip().split()
					type_ = int(line[2])
					my_data.atoms[type_]["q"] = float(line[-1])

	if args.dpkl is not None:
		with open(args.dpkl, 'wb') as f:
			pickle.dump(my_data, f)
	'''
	print(f"Atoms:")
	for key in my_data.atoms:
		typ = my_data.atoms[key]["type"]
		m   = my_data.atoms[key]["mass"]
		q   = my_data.atoms[key]["q"]
		eps = my_data.atoms[key]["eps"]
		sig = my_data.atoms[key]["sig"]
		print(f"atm_types[\"{typ}\"] = {{\"m\": {m}, \"q\": {q}, \"eps\": {eps}, \"sig\": {sig}}}")
	'''
	print(f"\n\nBonds:")
	for key in my_data.bonds:
		print(f"bond type: {key} with atom types: {my_data.bonds[key]['type']}")
		# for typ in my_data.bonds[key]["type"]:
			# kd = my_data.bonds[key]["kd"]
			# r0 = my_data.bonds[key]["r0"]
			# print(f"bon_types[{typ}] = {{\"type\": \"harmonic\", \"kb\": {kd}, \"r0\": {r0}}}")

	print(f"\n\nAngles:")
	for key in my_data.angles:
		print(f"angle type: {key} with atom types: {my_data.angles[key]['type']}")
		# for typ in my_data.angles[key]["type"]:
			# k0 = my_data.angles[key]["k0"]
			# t0 = my_data.angles[key]["t0"]
			# print(f"ang_types[{typ}] = {{\"type\": \"harmonic\", \"ka\": {k0}, \"a0\": {t0}}}")

	'''
	print(f"\n\nDihedrals:")
	for key in my_data.dih:
		for typ in my_data.dih[key]["type"]:
			V0 = my_data.dih[key]["V0"]
			V1 = my_data.dih[key]["V1"]
			V2 = my_data.dih[key]["V2"]
			V3 = my_data.dih[key]["V3"]
			print(f"dih_types[{typ}] = {{\"type\": \"opls\", \"kd\":[{V0}, {V1}, {V2}, {V3}]}}")

	print(f"\n\nImpropers:")
	for key in my_data.impr:
		for typ in my_data.impr[key]["type"]:
			k   = my_data.impr[key]["k"]
			phi = my_data.impr[key]["phi"]
			print(f"impr_types[{typ}] = {{\"type\": \"harmonic\", \"kb\": {k}, \"r0\": {phi}}}")
	'''
