#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import nglview as nv
import argparse
import pickle
import MDAnalysis.transformations as xforms
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(description="Compute hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--datafile",     dest="df",   type=str, action="store", help="enter address of datafile.",        required=True)
parser.add_argument("--settingsfile", dest="sf",   type=str, action="store", help="enter address of settingsfile.",    required=True)
parser.add_argument("--coords",       dest='c',    type=str, action="store", help="enter address of coordinate file.", required=True)
parser.add_argument("--pkl",          dest="pkl",  type=str, action="store", help="enter address of pickle dump.",     required=True)
parser.add_argument("--dpkl",         dest="dpkl",  type=str, action="store", help="enter address of pickle dump.",     required=True)
args = parser.parse_args()

# Define the possible neighbor displacements (including diagonals)
shifts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
							[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], 
							[1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1],
							[0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
							[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], 
							[-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]])

# get the lattice index 
def latticify(locs, li, grid_n):
	nx  = (locs[:,0] // li) % grid_n
	ny  = (locs[:,1] // li) % grid_n
	nz  = (locs[:,2] // li) % grid_n
	idx = np.array(nz * grid_n**2 + ny * grid_n + nx, dtype=int)
	return idx

# Helper function to get the 3D grid coordinates from a linear index
def linear_to_grid(j, grid_n):
	nz = j // (grid_n**2)
	ny = (j % (grid_n**2)) // grid_n
	nx = j % grid_n
	return np.array([nx, ny, nz])

# Helper function to convert grid coordinates back to a linear index
def grid_to_linear(n, grid_n):
	return int(n[2] * grid_n**2 + n[1] * grid_n + n[0])

# Function to get the neighbors of a grid point
def check_neighbors(occupied_cells, grid_n):

	check = True
	for j in occupied_cells:
		# Get the 3D coordinates of the grid point j
		coords = linear_to_grid(j, grid_n)

		# Loop through each shift and calculate the neighbor coordinates
		while check:
			for shift in shifts:
				# Apply the shift and apply modulo to handle periodic boundary conditions
				neighbor_coords = (coords + shift) % grid_n
				
				# Convert neighbor coordinates back to a linear index
				neighbor_idx = grid_to_linear(neighbor_coords, grid_n)

				# run the neighborhood check
				if neighbor_idx in occupied_cells:
					check = True
					break
				else:
					check = False
					continue

	return check

# get the appropriate cell_density 
def find_cell_density(universe, atom_types, grid_n=1):

	# get the atoms
	atom_selection = universe.select_atoms(atom_types)
	positions      = atom_selection.positions

	for loop in range(250):
		# define the cells array
		cells = np.zeros(grid_n**3)

		# get the subcube dimensions
		li = (np.prod(universe.dimensions[:3]) / (grid_n**3)) ** (1/3)

		# get the indices
		indices = latticify(positions, li, grid_n)

		# set up cells
		np.add.at(cells, indices, 1)

		# increase density of grids
		if np.any(cells > 8):
			grid_n += 1 # make finer  cells
		else:
			break

	# print(f"cells.shape = {cells.shape}, cells>1 = {cells>1}", flush=True)
	heavy_cells = np.array(range(grid_n**3),dtype=int)[cells > 1]
	light_cells = np.array(range(grid_n**3),dtype=int)[~(cells > 1)]
	return grid_n, heavy_cells, light_cells

if __name__=="__main__":

	# start the clock
	start = time.time()

	# get the names of the files
	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	# generate the universe object
	print(f"Making the universe...", end=' ', flush=True)
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"made!", flush=True)

	# unwrap the trajectory
	polymer = u.select_atoms("resid 1")
	other   = u.select_atoms("not resid 1")
	_all    = u.select_atoms("all")

	# set up the transforms
	transforms = [mda.transformations.unwrap(_all), mda.transformations.center_in_box(polymer, wrap=False), mda.transformations.wrap(_all)]
	u.trajectory.add_transformations(*transforms)

	# define polymer atom types 
	heavy_atoms = "type 1 or type 3-6 or type 8-9"

	# define the discretization object
	discretization = dict()

	# define a grid density
	grid_n = 1
	for idx, ts in enumerate(u.trajectory):
		print(f"ts = {ts.time} femtoseconds.", flush=True)
		atom_selection = u.select_atoms(heavy_atoms)
		grid_n, heavy_cells, light_cells = find_cell_density(u, heavy_atoms, grid_n)
		discretization[idx] = {"mesh_size": grid_n, "heavy_atoms": heavy_cells, "light_atoms": light_cells}
		grid_n = 1

	# create the grand object
	all_info = [u, discretization]

	file = open(args.pkl, 'wb')
	pickle.dump(all_info, file)
	file.close()

	"""
	print(f"Made all the objects!", flush=True)
	H_polymer = "type 7"
	H_water   = "type 11"
	O_polymer = "type 5"
	N_polymer = "type 6"
	O_water   = "type 10"

	# access atom types for hbonding
	print(f"Run hbond analysis...", flush=True)
	hbonds_with_Opoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=O_polymer, donors_sel=O_water)
	hbonds_with_Npoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=N_polymer, donors_sel=O_water)
	hbonds_with_Owater = HBA(universe=u, hydrogens_sel=H_polymer, acceptors_sel=O_water,   donors_sel=N_polymer)
	hbonds_with_Op_Np  = HBA(universe=u, hydrogens_sel=H_polymer, acceptors_sel=O_polymer, donors_sel=N_polymer)

	# run hbond analysis with MDA
	print(f"Running hbond analysis with water H being donated by water O and accepted by polymer O.", end=' ', flush=True)
	hbonds_with_Opoly.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with water H being donated by water O and accepted by polymer N.", end=' ', flush=True)
	hbonds_with_Npoly.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with amidic H being donated by amidic N and accepted by water O.", end=' ', flush=True)
	hbonds_with_Owater.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with amidic H being donated by amidic N and accepted by polymer O.", end=' ', flush=True)
	hbonds_with_Op_Np.run()
	print(f"done!", flush=True)
	"""

