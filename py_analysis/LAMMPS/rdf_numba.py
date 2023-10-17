#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re 
import argparse
import matplotlib.pyplot as plt
import time
import numba

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument("--trajfile", dest='traj', action='store', type=str, help="Name of trajfile.")
parser.add_argument("--atm1", dest="atm1", action='store', type=int, help="Name of atom type 1.")
parser.add_argument("--atm2", dest="atm2", action='store', type=int, help="Name of atom type 2.")
args = parser.parse_args()

v = [[0,0,0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = np.array([np.array(u) for u in v])

def create_traj_object(trajfile):

	traj = {}

	timestep_flag  = False
	numatoms_flag  = False
	boxbounds_flag = False
	itematoms_flag = False
	ts = None

	f = open(trajfile, 'r')
	for line in f:
		if re.findall(r'^ITEM: TIMESTEP$', line):
			timestep_flag = True
			numatoms_flag  = False
			boxbounds_flag = False
			itematoms_flag = False
			continue

		elif re.findall(r'^ITEM: NUMBER OF ATOMS$', line):
			timestep_flag  = False
			numatoms_flag  = True
			boxbounds_flag = False
			itematoms_flag = False
			continue

		elif re.findall(r'^ITEM: BOX BOUNDS pp pp pp$', line):
			timestep_flag  = False
			numatoms_flag  = False
			boxbounds_flag = True
			idx = 0
			itematoms_flag = False
			continue

		elif re.findall(r'^ITEM: ATOMS id type x y z', line):
			timestep_flag  = False
			numatoms_flag  = False
			boxbounds_flag = False
			itematoms_flag = True
			traj[ts]["coords"] = {}
			continue

		elif boxbounds_flag:
			contents = line.strip().split()
			if idx==0:
				traj[ts]["xlo"] = float(contents[0])
				traj[ts]["xhi"] = float(contents[1])
				idx+=1
			elif idx==1:
				traj[ts]["ylo"] = float(contents[0])
				traj[ts]["yhi"] = float(contents[1])
				idx+=1
			elif idx==2:
				traj[ts]["zlo"] = float(contents[0])
				traj[ts]["zhi"] = float(contents[1])
				idx+=1
			continue

		elif timestep_flag:
			if ts == None:
				pass
			else:
				print("Processed!", flush=True)
			contents = line.strip().split()
			ts = int(contents[0])
			traj[ts] = {}
			print(f"Processing timestep {ts}...", end=' ', flush=True)
			continue

		elif numatoms_flag:
			contents = line.strip().split()
			traj[ts]["natoms"] = contents[0]
			continue

		elif itematoms_flag:
			contents = line.strip().split()
			atm_sr   = int(contents[1])
			# figure out which molecule this particle is in
			if not(atm_sr in traj[ts]["coords"]):
				traj[ts]["coords"][atm_sr]   = np.empty((0,3))
			coords = np.array([float(contents[2]), float(contents[3]), float(contents[4])], dtype=np.float64)
			traj[ts]["coords"][atm_sr] = np.vstack((traj[ts]["coords"][atm_sr], coords))
			continue

	print ("Processed!", flush=True)
	n_atoms = {}
	for atm_sr in traj[ts]["coords"]:
		n_atoms[atm_sr] = len(traj[ts]["coords"][atm_sr])
	return traj, n_atoms

@numba.njit
def big_calcs(coords_1, coords_2, box_dims, rdf_total, atm_type_1, atm_type_2, N_1, N_2, num_bins):

	perturb_coords_2 = np.vstack([pdist+coords_2 for pdist in v])
	distances        = np.linalg.norm(coords_1[:,np.newaxis] - perturb_coords_2, axis=2)
	distances        = distances.reshape(-1)
	max_distance     = np.min(box_dims)/2.0
	distances        = distances[distances < max_distance]
	low_dists        = distances > 1e-2
	distances        = distances[low_dists]

	# create the histogram
	hist, bin_edges  = np.histogram(distances, bins=num_bins, range=(0, max_distance))

	# calculate the bin centers 
	bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
	bin_width   = bin_edges[1] - bin_edges[0]

	# density of particle 2
	rho_2 = N_2/np.prod(box_dims)
	hist = hist/(4*np.pi*(bin_centers**2)*bin_width*rho_2*N_1*27)

	rdf_total += hist

	return rdf_total

def calculate_rdf(traj, atm_type_1, atm_type_2, N_1, N_2, num_bins=1000):

	rdf_total = big_calcs(traj, atm_type_1, atm_type_2, N_1, N_2, num_bins)
	rdf_total /= len(traj)

	return rdf_total, bin_centers



if __name__=="__main__":

	start = time.time()
	# get all the information from the trajectory file
	traj, n_atoms = create_traj_object(args.traj)

	# now calculate the rdf
	rdf, bins = calculate_rdf(traj, args.atm1, args.atm2, n_atoms[args.atm1], n_atoms[args.atm2])

	fig = plt.figure(figsize=(5,5))
	ax  = plt.axes()

	ax.plot(bins, rdf, c='steelblue', lw=1)
	ax.set_xlim(0,10)
	ax.axhline(y=1, c='r', ls='--')
	fig.savefig('rdf', dpi=1200, bbox_inches="tight")
	stop = time.time()
	print(f"Time to calculate rdf is {stop-start} seconds.", flush=True)
