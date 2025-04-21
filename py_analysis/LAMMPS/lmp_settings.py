#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re

class Settings:
	def __init__(self, settingsfile):
		self.settingsfile = settingsfile
		self.pair_coeff  = dict()
		self.bond_coeff  = dict()
		self.angle_coeff = dict()
		self.dih_coeff   = dict()
		self.impr_coeff  = dict()
	
	def read_pairs(self, line):
		# self.pair_coeff[atom_type]  = (eps, sig)
		self.pair_coeff[int(line[1])] = (float(line[3]), float(line[4]))
		return

	def read_bonds(self, line):
		# self.bond_coeff[bond_type] = (k, r0)
		self.bond_coeff[int(line[1])] = (float(line[2]), float(line[3]))
		return

	def read_angles(self, line):
		# self.angle_coeff[angle_type] = (k, theta0)
		self.angle_coeff[int(line[1])] = (float(line[2]), float(line[3]))
		return

	def read_dih(self, line):
		# self.angle_coeff[dih_type] = (v0, v1, v2, v3)
		self.dih_coeff[int(line[1])] = (float(line[2]), float(line[3]), float(line[4]), float(line[5]))
		return 
	
	def read_impr(self, line):
		# self.impr_coeff[impr_type] = (k0, phi)
		self.impr_coeff[int(line[1])] = (float(line[2]), float(line[3]))
		return

	def read_settings(self):
		with open(self.settingsfile, 'r') as f:
			for line in f:
				if line.strip() == "":
					continue
				else:
					line_split = line.strip().split()
					if re.findall(r"^pair_coeff", line):
						self.read_pairs(line_split)
					elif re.findall(r"^bond_coeff", line):
						self.read_bonds(line_split)
					elif re.findall(r"^angle_coeff", line):
						self.read_angles(line_split)
					elif re.findall(r"^dihedral_coeff", line):
						self.read_dih(line_split)
					elif re.findall(r"^improper_coeff", line):
						self.read_impr(line_split)
		return
