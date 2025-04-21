#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re

class Data:
	def __init__(self, datafile):
		self.datafile     = datafile
		self.natoms       = 0
		self.nbonds       = 0
		self.nangles      = 0
		self.ndihedrals   = 0
		self.nimpropers   = 0
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
		self.read_data()
		return

	def check_section(self, line):
		if re.findall(r"^LAMMPS Description", line) or re.findall(r"^#LAMMPS", line):
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
		elif re.findall(r"^Velocities", line):
			self.b_description = False
			self.b_masses      = False
			self.b_atoms       = False
			self.b_bonds       = False
			self.b_angles      = False
			self.b_dihedrals   = False
			self.b_impropers   = False

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
		self.atoms[int(line[0])]["mass"]   = float(line[1])
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
		atoms_in_bond = tuple(sorted([self.atoms[atype_1]["type"], self.atoms[atype_2]["type"]]))
		self.bonds[btype]["type"].add(atoms_in_bond)
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
					elif self.b_angles or self.b_dihedrals or self.b_impropers:
						break
		return
