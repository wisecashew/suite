#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation

# Define the possible neighbor displacements (including diagonals)
shifts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
							[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0],
							[1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1],
							[0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
							[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
							[-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]], dtype=int)

norm_shifts = shifts / np.linalg.norm(shifts, axis=1)[:, np.newaxis]

# Helper function to get the 3D grid coordinates from a linear index
def linear_to_grid(j, grid_n):
	nz = j // (grid_n**2)
	ny = (j % (grid_n**2)) // grid_n
	nx = j % grid_n
	return np.array([nx, ny, nz])

# Helper function to convert grid coordinates back to a linear index
def grid_to_linear(n, grid_n):
	return int(n[2] * grid_n**2 + n[1] * grid_n + n[0])

# get lattice directions
def get_lat(mon2, mon1, occupied):
	diff           = mon2 - mon1
	projection     = np.dot(norm_shifts, diff)
	sorted_indices = np.argsort(projection)[::-1]
	counter        = 0
	for i in sorted_indices:
		best_dir = shifts[i]
		exists = np.any(np.all(np.array(occupied) == occupied[-1]+best_dir, axis=1))
		if exists:
			counter += 1
			continue
		else:
			occupied.append(occupied[-1]+best_dir)
			break
	if counter == 26:
		print(f"Not possible to latticify!")
		exit()
	return 

def plot_cube(ax, center, size=1, color="blue"):
	x, y, z = center
	half    = size / 2

	# Vertices of the cube
	vertices = np.array([
		[x - half, y - half, z - half],
		[x + half, y - half, z - half],
		[x + half, y + half, z - half],
		[x - half, y + half, z - half],
		[x - half, y - half, z + half],
		[x + half, y - half, z + half],
		[x + half, y + half, z + half],
		[x - half, y + half, z + half]
	])

	# Define the six faces of the cube
	faces = [
		[vertices[j] for j in [0, 1, 2, 3]],  # Bottom
		[vertices[j] for j in [4, 5, 6, 7]],  # Top
		[vertices[j] for j in [0, 3, 7, 4]],  # Left
		[vertices[j] for j in [1, 2, 6, 5]],  # Right
		[vertices[j] for j in [0, 1, 5, 4]],  # Front
		[vertices[j] for j in [3, 2, 6, 7]]   # Back
	]

	# Add the faces to the plot
	ax.add_collection3d(Poly3DCollection(faces, facecolors=color, edgecolors='black', linewidths=0.5, alpha=0.8, zorder=1))
	return

def plot_bonds(ax, positions):
	"""
	Draws bonds (lines) between consecutive cubes based on their positions.
	"""
	for i in range(len(positions) - 1):
		x = [positions[i][0] + 0.5, positions[i + 1][0] + 0.5]  # Midpoints of cubes
		y = [positions[i][1] + 0.5, positions[i + 1][1] + 0.5]
		z = [positions[i][2] + 0.5, positions[i + 1][2] + 0.5]
		ax.plot(x, y, z, color="black", linewidth=2)
	return

class Lattice_LMP:
	def __init__(self, universe, hba):
		self.universe = universe
		self.hba      = hba
		self.lattice  = dict()
		return

	# get the start and end indices
	get_start = lambda self, m: 1  if m == 0 else 21 + 19 * (m-1)
	get_end   = lambda self, m: self.get_start(m)+19 if (m == 0 or m == 39) else self.get_start(m)+18

	def add_monomer_to_lattice(self, idx, m_to_add, m_previous):
		diff           = m_to_add - m_previous
		projection     = np.dot(norm_shifts, diff)
		sorted_indices = np.argsort(projection)[::-1]
		counter        = 0
		for i in sorted_indices:
			best_dir = shifts[i]
			exists = np.any(np.all(np.array(self.lattice[idx]["all"]) == self.lattice[idx]["polymer"][-1]+best_dir, axis=1))
			if exists:
				print(f"{counter}: Encountered particle while adding monomer. pushing off...", flush=True)
				counter += 1
				continue
			else:
				self.lattice[idx]["polymer"] = np.vstack((self.lattice[idx]["polymer"], self.lattice[idx]["polymer"][-1]+best_dir))
				self.lattice[idx]["all"]     = np.vstack((self.lattice[idx]["all"],     self.lattice[idx]["polymer"][-1]+best_dir))
				break
		if counter == 26:
			print(f"Not possible to latticify while adding monomers!", flush=True)
			# exit()
		return 

	def add_solvent_to_lattice(self, idx, s_to_add, m_relevant):
		diff           = s_to_add - m_relevant
		projection     = np.dot(norm_shifts, diff)
		sorted_indices = np.argsort(projection)[::-1]
		counter        = 0
		for i in sorted_indices:
			best_dir = shifts[i]
			exists = np.any(np.all(np.array(self.lattice[idx]["all"]) == m_relevant+best_dir, axis=1))
			if exists:
				print(f"{counter}: Encountered particle while adding solvent. pushing off...", flush=True)
				counter += 1
				continue
			else:
				self.lattice[idx]["solvent"] = np.vstack((self.lattice[idx]["solvent"], m_relevant+best_dir))
				self.lattice[idx]["all"]     = np.vstack((self.lattice[idx]["all"],     m_relevant+best_dir))
				break
		if counter == 26:
			print(f"Not possible to add while adding solvents!", flush=True)
			# exit()
		return 

	def latticify_dry(self):
		c_of_m = list()
		for idx,ts in enumerate(self.universe.trajectory):
			print(f"At timestep {ts.time}.", flush=True)
			c_of_m.clear()
			self.lattice[idx] = dict()
			self.lattice[idx]["all"]     = np.empty((0,3), dtype=int)
			self.lattice[idx]["polymer"] = np.empty((0,3), dtype=int)
			for midx in range(40):
				sel_str    = f"bynum {self.get_start(midx)}:{self.get_end(midx)}"
				selection  = self.universe.select_atoms(sel_str)
				c_of_m.append(selection.center_of_mass())
				if midx == 0:
					self.lattice[idx]["polymer"] = np.vstack((self.lattice[idx]["polymer"], np.array([0,0,0], dtype=int)))
				else:
					self.add_monomer_to_lattice(idx, c_of_m[-1], c_of_m[-2])
			
		return

	def latticify(self):
		c_of_m = list()
		w_mols = list()
		for idx,ts in enumerate(self.universe.trajectory):
			print(f"At timestep {ts.time}.", flush=True)
			c_of_m.clear()
			w_mols.clear()
			self.lattice[idx] = dict()
			self.lattice[idx]["all"]     = np.empty((0,3), dtype=int)
			self.lattice[idx]["polymer"] = np.empty((0,3), dtype=int)
			self.lattice[idx]["solvent"] = np.empty((0,3), dtype=int)
			for midx in range(40):
				sel_str    = f"bynum {self.get_start(midx)}:{self.get_end(midx)}"
				selection  = self.universe.select_atoms(sel_str)
				c_of_m.append(selection.center_of_mass())
				if midx == 0:
					self.lattice[idx]["polymer"] = np.vstack((self.lattice[idx]["polymer"], np.array([0,0,0], dtype=int)))
				else:
					self.add_monomer_to_lattice(idx, c_of_m[-1], c_of_m[-2])

			# update lattice molecules
			self.lattice[idx]["all"] = np.vstack((self.lattice[idx]["all"], self.lattice[idx]["polymer"]))
			
			# get the hydrogen bonds at this timestep
			hb_Ow  = self.hba["Op_with_water"][0][self.hba["Op_with_water"][0][:,0] == idx]
			hb_Nw  = self.hba["Np_with_water"][0][self.hba["Np_with_water"][0][:,0] == idx]
			hb_net = np.vstack((hb_Ow, hb_Nw))

			# water oxygens participating in h-bonding
			water_oxygen_idx = hb_net[:, 1]
			_, indices = np.unique(water_oxygen_idx, return_index=True)
			hb_net = hb_net[np.sort(indices)]

			# get the water molecules that are H-bonding
			for hb in hb_net:
				midx = (hb[-1]-2)//19 if (hb[-1]-2)//19 < 40 else 39
				water = self.universe.select_atoms(f"bynum {hb[1]} or bynum {hb[2]}")
				resid = water.residues.resids[0]
				w_mols.append((midx, resid))

			# get the water molecules along with the resid
			w_mols = list(set(w_mols))

			# loop through each water molecule
			for w in w_mols:
				water = self.universe.select_atoms(f"resid {w[1]}")
				self.add_solvent_to_lattice(idx, water.center_of_mass(), self.lattice[idx]["polymer"][w[0]])
			
		return

	def latticify_v2(self):
		c_of_m = list()
		w_mols = list()
		for idx, ts in enumerate(self.universe.trajectory):
			print(f"At timestep {ts.time}.", flush=True)
			c_of_m.clear()
			w_mols.clear()
			self.lattice[idx] = dict()
			self.lattice[idx]["all"]     = np.empty((0,3), dtype=int)
			self.lattice[idx]["polymer"] = np.empty((0,3), dtype=int)
			self.lattice[idx]["solvent"] = np.empty((0,3), dtype=int)

			# get the hydrogen bonding info and get rid of duplicates
			hb_Ow  = self.hba["Op_with_water"][0][self.hba["Op_with_water"][0][:,0] == idx]
			hb_Nw  = self.hba["Np_with_water"][0][self.hba["Np_with_water"][0][:,0] == idx]
			hb_net = np.vstack((hb_Ow, hb_Nw))

			# water oxygens participating in h-bonding
			water_oxygen_idx = hb_net[:, 1]
			_, indices = np.unique(water_oxygen_idx, return_index=True)
			hb_net = hb_net[np.sort(indices)]
			# print(hb_net)

			# get the water molecules that are H-bonding
			w_mols = [[] for _ in range(40)]
			for hb in hb_net:
				midx = (hb[-1]-2)//19 if (hb[-1]-2)//19 < 40 else 39
				water = self.universe.select_atoms(f"bynum {hb[1]} or bynum {hb[2]}")
				resid = water.residues.resids[0]
				w_mols[midx].append(resid)
				# print(w_mols)
			# print(w_mols)
			# exit()
			for midx in range(40):
				sel_str    = f"bynum {self.get_start(midx)}:{self.get_end(midx)}"
				selection  = self.universe.select_atoms(sel_str)
				c_of_m.append(selection.center_of_mass())
				# add the monomer to the lattice
				if midx == 0:
					self.lattice[idx]["polymer"] = np.vstack((self.lattice[idx]["polymer"], np.array([0,0,0], dtype=int)))
				else:
					self.add_monomer_to_lattice(idx, c_of_m[-1], c_of_m[-2])
				# add the solvents to the lattice
				for mol in w_mols[midx]:
					water = self.universe.select_atoms(f"resid {mol}")
					self.add_solvent_to_lattice(idx, water.center_of_mass(), self.lattice[idx]["polymer"][midx])
			
		return
	
	def get_Nmm(self, ts):

		# define some variables 
		Nmm     = -1
		polymer = self.lattice[ts]["polymer"]

		# compute pairwise differences
		diff = polymer[:, np.newaxis, :] - polymer[np.newaxis, :, :]

		# flatten differences to check against shifts
		diff = diff.reshape(-1, 3)

		# check with shifts
		matches = np.any(np.all(diff[:, np.newaxis, :] == shifts, axis=2), axis=1)
		
		# count matches
		Nmm = np.sum(matches) // 2

		# update
		self.lattice[ts]["Nmm"] = Nmm

		return
	
	def get_Nmsa(self, ts):

		# define some variables
		Nmsa = -1
		polymer  = self.lattice[ts]["polymer"]
		solvents = self.lattice[ts]["solvent"]

		# compute pairwise differences
		diff = polymer[:, np.newaxis, :] - solvents[np.newaxis, :, :]
		diff = diff.reshape(-1, 3)

		# check with shifts
		matches = np.any(np.all(diff[:, np.newaxis, :] == shifts, axis=2), axis=1)

		# count matches 
		Nmsa = int(np.sum(matches))

		# update
		self.lattice[ts]["Nmsa"] = Nmsa

		return

