#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import argparse
import pickle 
import argparse

parser = argparse.ArgumentParser(description="Add the impropers to a pnipam datafile.")
parser.add_argument("--datafile",      dest='df',    type=str, action='store', help="enter address of datafile.")
parser.add_argument("--ndatafile",     dest='ndf',   type=str, action='store', help="enter address of new datafile.")
parser.add_argument("--settingsfile",  dest='sf',    type=str, action='store', help="enter address of settingsfile.")
parser.add_argument("--imp_sr",        dest='impsr', type=int, action='store', help="enter indices of the important atom types.", nargs='+')
args = parser.parse_args()

class Data:
	def __init__(self, datafile):
		self.datafile   = datafile
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
			self.bonds[btype]["atoms"]   = list()
		self.bonds[btype]["atoms"].append((int(line[2]), int(line[3])))
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

	def read_data(self):
		with open(self.datafile, 'r') as f:
			for line in f:
				if line.strip() == "":
					continue
				else:
					self.check_section(line)
					line_split = line.strip().split()
					if re.findall(r"[a-zA-Z]+", line_split[0]):
						continue
					elif self.b_description:
						self.read_description(line)
					elif self.b_masses:
						self.read_masses(line)
					elif self.b_atoms:
						self.read_atoms(line)
					elif self.b_bonds:
						self.read_bonds(line)
					elif self.b_angles:
						break
		return

if __name__=="__main__":
	
	my_data = Data(args.df)

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
					my_data.read_bonds(line)
				elif my_data.b_angles:
					break


	impr_5 = set()
	impr_7 = set()

	# get all the particle serial numbers of type 5
	n5_list = list()
	n7_list = list()

	for sr_no in my_data.atoms_serial:
		if my_data.atoms_serial[sr_no] == args.impsr[0]:
			n5_list.append(sr_no)
		elif my_data.atoms_serial[sr_no] == args.impsr[1]:
			n7_list.append(sr_no)
	
	for n5 in n5_list:
		impr = [n5]
		for btype in my_data.bonds:
			for atoms in my_data.bonds[btype]["atoms"]:
				if n5 == atoms[0]:
					impr.append(atoms[1])
				elif n5 == atoms[1]:
					impr.append(atoms[0])
				else:
					continue
		impr_5.add(tuple(impr))
		if len(impr) == 4:
			pass
		else:
			print(f"Something is fucked.")

	for n7 in n7_list:
		impr = [n7]
		for btype in my_data.bonds:
			for atoms in my_data.bonds[btype]["atoms"]:
				if n7 == atoms[0]:
					impr.append(atoms[1])
				elif n7 == atoms[1]:
					impr.append(atoms[0])
				else:
					continue
		impr_7.add(tuple(impr))
		if len(impr) == 4:
			pass
		else:
			print(f"Something is fucked.")

	# write the improper
	with open(args.df, 'a') as f:
		f.write("\n")
		f.write("Impropers\n")
		f.write("\n")
		idx = 1
		for impr in impr_5:
			f.write(f"{idx} 1 {impr[0]} {impr[1]} {impr[2]} {impr[3]}\n")
			idx += 1
		
		for impr in impr_7:
			f.write(f"{idx} 2 {impr[0]} {impr[1]} {impr[2]} {impr[3]}\n")
			idx += 1


	f = open(args.ndf, 'w')
	g = open(args.df, 'r')

	for line in g:
		f.write(line)
		if re.findall(r"dihedral types", line):
			f.write(f"{len(impr_5) + len(impr_7)} impropers\n")
			f.write(f"2 improper types\n")

	f.close()
	g.close()

	g = open(args.sf, 'a')
	g.write("########################################\n")
	g.write("########################################\n")
	g.write("improper_coeff 1 10.5 180\n")
	g.write("improper_coeff 2 1.0  180")
