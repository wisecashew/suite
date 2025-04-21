#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import argparse

parser = argparse.ArgumentParser(description="Process hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--pkl",          dest="pkl",   type=str, action="store", required=False, help="enter address of pickle file to store H-bonding info.")
args = parser.parse_args()

# Define the possible neighbor displacements (including diagonals)
shifts = np.array([[0,0,0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
							[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], 
							[1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1],
							[0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
							[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], 
							[-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]], dtype=int)

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

def check(latt):
	if (latt[0] - latt[1] in shifts) and (latt[0] - latt[2] in shifts) and (latt[1] - latt[2] in shifts):
		return True
	else:
		return False

if __name__=="__main__":

	# set up the start
	start = time.time()

	# get the h-bonding info
	hb_file = open(f"TOT_HBONDS/hbond.0.pkl", "rb")
	hbs     = pickle.load(hb_file)
	hb_file.close()

	# get the centering info
	traj_file = open(f"TOT_CLEANED_COORDS/CENTERED/universe.0.pkl", "rb")
	traj      = pickle.load(traj_file)
	traj_file.close()

	# get the universe and discretization
	u              = traj[0]
	discretization = traj[1]

	ATOMIC_INFO = list()
	GEOM_INFO   = list()
	# LATT_INFO   = dict()

	# access the results
	for i in range(4):
		ATOMIC_INFO.append(np.array(hbs[i].results.hbonds[:,0:4], dtype=int))
		GEOM_INFO.append(np.array(hbs[i].results.hbonds[:,4:], dtype=int))
		# LATT_INFO[i] = dict()

	for idx, ts in enumerate(u.trajectory):
		print(f"At timestep = {ts.time}.", flush=True)

		# get the all locations
		locations = u.select_atoms("all").positions

		# get the subcube length
		mesh_size = discretization[idx]["mesh_size"]
		li = (np.prod(u.dimensions[:3]) / (mesh_size**3)) ** (1/3)

		# get the grid point for the donor oxygen, water hydrogen, and polymer oxygen
		for i in range(4):
			for a_info in ATOMIC_INFO[i][ATOMIC_INFO[i][:,0] == idx]:
				atomic_locs = np.array([locations[a_info[1]], locations[a_info[2]], locations[a_info[3]]])
				latts       = latticify(atomic_locs, li, discretization[idx]["mesh_size"])
				grid        = np.array([linear_to_grid(latts[0], mesh_size), linear_to_grid(latts[1], mesh_size), linear_to_grid(latts[2], mesh_size)], dtype=int)
				if not check(grid):
					print(f"Something is too distant @ {i}!")
					print(grid)
					exit()
				# LATT_INFO[i][idx].append((latts, grid))

	# set the stop
	stop = time.time()
	print(f"Computation took {stop-start} seconds.", flush=True)

