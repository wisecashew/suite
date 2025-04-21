#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re

class Configurations:
	def __init__(self, coordsfile, DELTA_TS):
		self.coordsfile  = coordsfile
		self.delta_ts    = DELTA_TS
		self.timesteps   = []
		self.natoms      = []
		self.box_dims    = []
		self.xdim        = []
		self.ydim        = []
		self.zdim        = []
		self.coords      = []
		self.r_timestep  = r"ITEM: TIMESTEP"
		self.timestep_b  = False
		self.r_natoms    = r"ITEM: NUMBER OF ATOMS"
		self.natoms_b    = False
		self.r_boxdims   = r"ITEM: BOX BOUNDS pp pp pp"
		self.boxdims_b   = False
		self.r_colnames  = r"ITEM: ATOMS id type x y z"
		self.cols_b      = False
		self.consider_ts = True
		self.read_coordsfile()
		return

	def read_coordsfile(self):
		cfile = open(self.coordsfile, 'r')
		for line in cfile:
			if re.match(self.r_timestep, line):
				self.cols_b     = False 
				self.timestep_b = True
				try:
					if len(coords) > 0:
						self.coords.append(coords)
						coords = []
					else:
						pass
				except NameError:
					pass
			elif self.timestep_b:
				current_timestep = int(line.strip().split()[0])
				if len(self.timesteps) == 0:
					self.timesteps.append(current_timestep)
					self.timestep_b  = False
					self.consider_ts = True
				elif current_timestep == 0:
					self.timestep_b  = False
					self.consider_ts = False
				elif current_timestep < self.timesteps[-1]:
					self.timesteps.append(self.timesteps[-1] + self.delta_ts)
					self.timestep_b  = False
					self.consider_ts = True
				elif current_timestep > self.timesteps[-1]:
					self.timesteps.append(current_timestep)
					self.timestep_b  = False
					self.consider_ts = True
				print(f"At timestep {self.timesteps[-1]}.", flush=True)
			elif self.consider_ts:
				if re.match(self.r_natoms, line):
					self.natoms_b = True 
				elif self.natoms_b:
					self.natoms.append(int(line.strip().split()[0]))
					self.natoms_b = False
				elif re.match(self.r_boxdims, line):
					self.boxdims_b = True 
					dim_idx = 0
				elif self.boxdims_b:
					properties = line.strip().split()
					if dim_idx == 0:
						self.xdim.append((float(properties[0]), float(properties[1])))
						dim_idx += 1
					elif dim_idx == 1:
						self.ydim.append((float(properties[0]), float(properties[1])))
						dim_idx += 1
					elif dim_idx == 2:
						self.zdim.append((float(properties[0]), float(properties[1])))
						self.boxdims_b = False
				elif re.match(self.r_colnames, line):
					self.cols_b = True
					coords = []
				elif self.cols_b:
					properties = line.strip().split()
					coords.append([int(properties[0]), int(properties[1]), float(properties[2]), float(properties[3]), float(properties[4])])
		cfile.close()

		self.coords.append(coords.copy())
		coords.clear() 

		return
	
