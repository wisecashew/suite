import numpy as np
import re 
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument("--trajfile", dest='traj', action='store', type=str, help="Name of trajfile.")
parser.add_argument("--atm1", dest="atm1", action='store', type=int, help="Name of atom type 1.")
parser.add_argument("--atm2", dest="atm2", action='store', type=int, help="Name of atom type 2.")
args = parser.parse_args()

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
			print(f"ts = {ts}.")
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
	n_atoms = []
	for atm_sr in traj[ts]["coords"]:
		n_atoms.append( len(traj[ts]["coords"][atm_sr]) )
	return traj, n_atoms

def calculate_rdf(traj, atm_type_1, atm_type_2, N_1, N_2, num_bins=1000):

	"""
	Calculate RDF between particles of type 'a' and 'b' considering periodic boundary conditions using vectorization.

	Parameters:
	- coordinates_a: Numpy array of coordinates of particles of type 'a'.
	- coordinates_b: Numpy array of coordinates of particles of type 'b'.
	- atm_type_1: Particle type for 1.
	- atm_type_2: Particle type for 2.

	Returns:
	- rdf: RDF values for each bin.
	- bin_centers: Radial distance bin centers.
	"""
	rdf_total = np.zeros(num_bins, dtype=np.float64)

	# Pairwise distances (using broadcasting) between particles of type 'a' and 'b'
	for ts in traj:
		box_dims = np.array([ traj[ts]["xhi"]-traj[ts]["xlo"], traj[ts]["yhi"]-traj[ts]["ylo"], traj[ts]["zhi"]-traj[ts]["zlo"] ], dtype=np.float64)
		coords_1 = traj[ts]["coords"][atm_type_1]
		coords_2 = traj[ts]["coords"][atm_type_2]
		delta_r  = coords_1[:, np.newaxis, :] - coords_2
		delta_r -= box_dims * np.round(delta_r / box_dims)

		r = np.sqrt(np.sum(delta_r**2, axis=2))

		max_distance = np.min(box_dims)/2.0
		bin_edges = np.linspace(0, max_distance, num_bins + 1)

		# calculate rdf
		rdf, _ = np.histogram(delta_r, bins=bin_edges)
		normalization_factor = N_1 / (4*np.pi*bin_edges**2*(bin_edges[1]-bin_edges[0]))
		rdf = np.array(rdf, dtype=np.float64)
		rdf /= normalization_factor
		bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
		rdf_total += rdf

	rdf_total /= len(traj)

	return rdf_total, bin_centers



if __name__=="__main__":

	# get all the information from the trajectory file
	traj, n_atoms = create_traj_object(args.traj)
	print(traj)

	# now calculate the rdf
	rdf, bins = calculate_rdf(traj, args.atm1, args.atm2, n_atoms[0], n_atoms[1])

	fig = plt.figure(figsize=(5,5))
	ax  = plt.axes()

	ax.plot(bins, rdf, c='coral', marker='o', mec='k')
	fig.savefig('rdf', dpi=1200, bbox_inches="tight")

