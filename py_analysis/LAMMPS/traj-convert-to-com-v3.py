#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np 
import re
import argparse
import time
import copy
import pickle
import multiprocessing
multiprocessing.set_start_method('fork')
import itertools
import sys

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument ("--trajfile", dest='traj', action='store', type=str, help="Name of trajfile.")
parser.add_argument ("--ntrajfile", dest='ntraj', action='store', type=str, help="Name of cg'd trajfile.")
parser.add_argument ("--topology", dest="top", action="store", type=str, help="Enter name of pickled simulation topology.", default=None)
parser.add_argument ("--trajectory", dest="trajpkl", action="store", type=str, help="Enter name of pickled trajectory dictionary.", default=None)
parser.add_argument ("--nproc", dest='nproc', action='store', type=int, help="Number of processes to request.")
parser.add_argument ("--datafile", dest='data', action='store', type=str, help="Name of original datafile.")

args = parser.parse_args()

def create_simulation_object(datafile):

	# read the data file 
	sim_data = dict()

	extract_masses = False
	extract_atoms  = False
	extract_bonds  = False


	f = open(datafile, 'r')

	for line in f:
		if re.findall(r'atoms$', line):
			sim_data["natoms"] = int((line.strip().split())[0])
			continue
		elif re.findall(r'atom types$', line): 
			sim_data["natom_types"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'bonds$', line):
			sim_data["nbonds"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'bond types$', line): 
			sim_data["nbond_types"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'angles$', line):
			sim_data["nangles"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'angle types$', line): 
			sim_data["nangle_types"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'dihedrals$', line):
			sim_data["nangles"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'dihedral types$', line): 
			sim_data["nangle_types"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'impropers$', line):
			sim_data["nangles"] = int((line.strip().split()[0]))
			continue
		elif re.findall(r'improper types$', line): 
			sim_data["nangle_types"] = int((line.strip().split()[0]))
			continue
		elif not line.strip().split():
			continue
		elif re.findall(r'^Masses$', line):
			extract_atoms  = False
			extract_bonds  = False
			extract_masses = True
			sim_data["masses"] = {}
			continue
		elif re.findall(r'^Atoms$', line):
			extract_atoms  = True
			extract_bonds  = False
			extract_masses = False
			sim_data["atoms"] = {}
			sim_data["molecules"] = {}
			continue
		elif re.findall(r'^Bonds$', line):
			extract_atoms  = False
			extract_bonds  = True
			extract_masses = False
			sim_data["bonds"] = []
			continue

		elif re.findall(r'^Angles$', line):
			extract_atoms  = False
			extract_bonds  = False
			extract_masses = False
			continue

		elif re.findall(r'^Dihedrals$', line):
			extract_atoms  = False
			extract_bonds  = False
			extract_masses = False
			continue

		elif re.findall(r'^Impropers$', line):
			extract_atoms  = False
			extract_bonds  = False
			extract_masses = False
			continue

		elif extract_masses:
			content = line.strip().split()
			sim_data["masses"][int(content[0])] = float(content[1])
			continue

		elif extract_atoms:
			content = line.strip().split()
			sim_data["atoms"][int(content[0])] = {}
			sim_data["atoms"][int(content[0])]["molid"]   = int(content[1])
			sim_data["atoms"][int(content[0])]["atm_num"] = int(content[2])
			sim_data["atoms"][int(content[0])]["q"]       = float(content[3])
			if not(int(content[1]) in sim_data["molecules"]):
				sim_data["molecules"][int(content[1])] = {}
				sim_data["molecules"][int(content[1])]["atoms"]     = []
				sim_data["molecules"][int(content[1])]["atom_type"] = []
				sim_data["molecules"][int(content[1])]["bonds"]     = {}
			sim_data["molecules"][int(content[1])]["atoms"].append(int(content[0]))
			sim_data["molecules"][int(content[1])]["atom_type"].append(int(content[2]))
			sim_data["molecules"][int(content[1])]["bonds"][int(content[0])] = []
			continue

		elif extract_bonds:
			content = line.strip().split()
			sim_data["bonds"].append((int(content[2]), int(content[3])))
			# search which molecule the atoms are in 
			for molid in sim_data["molecules"]:
				if (int(content[2]) in sim_data["molecules"][molid]["atoms"]):
					sim_data["molecules"][molid]["bonds"][int(content[2])].append(int(content[3]))
					sim_data["molecules"][molid]["bonds"][int(content[3])].append(int(content[2]))
					break
	sim_data["nmolecules"] = len(sim_data["molecules"])

	# make tuple out of the list of atom types in sim_data["molecules"][molid]["atom_type"]
	for molid in sim_data["molecules"]:
		sim_data["molecules"][molid]["atom_type"] = tuple(sim_data["molecules"][molid]["atom_type"])

	# find number of different types of molecules
	mol_set  = set()
	ATM_LIST = []
	mol_tag  = 0
	for srno in sim_data["atoms"]:
		if mol_tag == 0:
			mol_tag = sim_data["atoms"][srno]["molid"]
		elif sim_data["atoms"][srno]["molid"] == mol_tag:
			pass
		else:
			# update mol_set with the previous ATM_LIST
			mol_set.add(tuple(ATM_LIST))
			ATM_LIST.clear()
			# clear the list so we can start plugging in the atoms for the next molecule
			# define the new mol_tag
			mol_tag = sim_data["atoms"][srno]["molid"]
		ATM_LIST.append(sim_data["atoms"][srno]["atm_num"])
	sim_data["nspecies"] = len(mol_set)

	# now that i know molecules exist, figure out which molid is which molecule
	mol_name = dict()
	for idx,element in enumerate(mol_set):
		mol_name[idx+1] = element

	for molid in sim_data["molecules"]:
		for molname in mol_name:
			if sim_data["molecules"][molid]["atom_type"] == mol_name[molname]:
				sim_data["molecules"][molid]["mol_name"] = molname
				break
			else:
				continue

	return sim_data



def divide_list_into_n_slices(list_to_divide, n_slices):
	"""Divides a list into N slices.

	Args:
	list_to_divide: The list to divide.
	n_slices: The number of slices to divide the list into.

	Returns:
	A list of lists, where each sublist contains the elements in one slice of the
	original list.
	"""

	# Check if the list is empty.
	if len(list_to_divide) == 0:
		return []

	# Check if the number of slices is greater than the length of the list.
	if n_slices > len(list_to_divide):
		raise ValueError("Number of slices must be less than or equal to the length of the list.")

	# Calculate the slice size.
	slice_size = len(list_to_divide) // n_slices

	# Create a list to store the slices.
	slices = []

	# Iterate over the slices and add the elements in each slice to the slices list.
	for i in range(n_slices):
		start_index = i * slice_size
		if i == n_slices - 1:
			sl = list_to_divide[start_index:]
		else:
			end_index = (i + 1) * slice_size
			sl = list_to_divide[start_index:end_index]
		slices.append(sl)

	# Return the list of slices.
	return slices

def traj_poker(trajfile):
	timesteps = []
	timestep_flag  = False
	ts = None

	f = open(trajfile, 'r')
	for line in f:
		if re.findall(r'^ITEM: TIMESTEP$', line):
			timestep_flag = True
			continue
		elif timestep_flag:
			contents = line.strip().split()
			timesteps.append(int(contents[0]))
			timestep_flag = False
			continue
	
	return timesteps



def traj_prober(sim_info, trajfile, timestep):
	traj = dict()
	approved       = False
	timestep_flag  = False
	numatoms_flag  = False
	boxbounds_flag = False
	itematoms_flag = False
	ts = None

	f = open(trajfile, 'r')
	for line in f:
		if re.findall(r'^ITEM: TIMESTEP$', line):
			timestep_flag  = True
			numatoms_flag  = False
			boxbounds_flag = False
			itematoms_flag = False
			continue

		elif approved and re.findall(r'^ITEM: NUMBER OF ATOMS$', line):
			timestep_flag  = False
			numatoms_flag  = True
			boxbounds_flag = False
			itematoms_flag = False
			continue

		elif approved and re.findall(r'^ITEM: BOX BOUNDS pp pp pp$', line):
			timestep_flag  = False
			numatoms_flag  = False
			boxbounds_flag = True
			idx = 0
			itematoms_flag = False
			continue

		elif approved and re.findall(r'^ITEM: ATOMS id type x y z', line):
			timestep_flag  = False
			numatoms_flag  = False
			boxbounds_flag = False
			itematoms_flag = True
			continue

		elif approved and boxbounds_flag:
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
			print(line)
			contents = line.strip().split()
			ts = int(contents[0])
			if ts in timestep:
				approved = True
				contents = line.strip().split()
				traj[ts] = {}
				print(f"Processing timestep {ts}...", end=' ', flush=True)
				continue
			else:
				approved = False
				continue

		elif numatoms_flag:
			contents = line.strip().split()
			traj[ts]["natoms"] = contents[0]
			traj[ts]["molid"]  = {}
			continue

		elif itematoms_flag:
			contents = line.strip().split()
			atm_sr   = int(contents[0])
			# figure out which molecule this particle is in
			for molid in sim_info["molecules"]:
				if atm_sr in sim_info["molecules"][molid]["atoms"]:
					if not(molid in traj[ts]["molid"]):
						traj[ts]["molid"][molid] = {}
						traj[ts]["molid"][molid]["coords"]   = np.empty((0,3))
						traj[ts]["molid"][molid]["bond_map"] = sim_info["molecules"][molid]["bonds"]
						traj[ts]["molid"][molid]["masses"]   = np.empty(0, dtype=np.float64)
						traj[ts]["molid"][molid]["mol_name"] = sim_info["molecules"][molid]["mol_name"]

					coords = np.array([float(contents[2]), float(contents[3]), float(contents[4])], dtype=np.float64)
					traj[ts]["molid"][molid]["coords"] = np.vstack((traj[ts]["molid"][molid]["coords"], coords))
					traj[ts]["molid"][molid]["masses"] = np.hstack((traj[ts]["molid"][molid]["masses"], sim_info["masses"][int(contents[1])]))
					break
				else:
					continue
			continue
	f.close()
	return traj


def create_traj_object(sim_info, trajfile, timesteps, nproc):

	# traj_prober has been defined above
	# time to send out the processes

	ts_list = divide_list_into_n_slices(timesteps, nproc)
	pool = multiprocessing.Pool(processes=nproc)
	results = pool.starmap(traj_prober, zip(itertools.repeat(sim_info), itertools.repeat(trajfile), ts_list))
	pool.close()
	pool.join()

	print(results)

	traj = results[0]
	for idx,d in enumerate(results):
		if idx == 0:
			continue
		else:
			traj.update(d)

	return traj

def unwrap_molecule(molecule, bonds, box_dims, sr_llim):

	untested = list(range(sr_llim, sr_llim+len(molecule)))
	tested   = []
	queue    = []

	while untested:
		wait = []
		if not queue:
			queue.append(untested[0])
		for i in queue:
			# print(f"i = {i}")
			neighbors = bonds[i]
			neighbors = [ni for ni in neighbors if ni not in tested]
			ri = molecule[i-sr_llim]
			for j in neighbors:
				rj = molecule[j-sr_llim]
				dr = rj - ri
				shift = np.round(dr/box_dims)
				molecule[j-sr_llim] -= shift*box_dims
				# print(f"j = {j}")
				# print(f"bonds[j] = {bonds[j]}")
				bonds[j].remove(i)
			tested.append(i)
			# print(f"right before remove, i = {i}...")
			untested.remove(i)
			wait.extend(neighbors)
		# print(f"queue = {queue}")
		queue = list(set(wait[:]))

	return # molecule

def unwrap_trajectory(traj_info):

	for ts in traj_info:
		print(f"Unwrapping timestep = {ts}...", end=' ', flush=True)
		sr_llim = 1
		for molid in traj_info[ts]["molid"]:
			box_dims = np.array([traj_info[ts]["xhi"]-traj_info[ts]["xlo"], \
				traj_info[ts]["yhi"]-traj_info[ts]["ylo"], \
				traj_info[ts]["zhi"]-traj_info[ts]["zlo"]])
			traj_info[ts]["molid"][molid]["coords"][:,0] -= traj_info[ts]["xlo"]
			traj_info[ts]["molid"][molid]["coords"][:,1] -= traj_info[ts]["ylo"]
			traj_info[ts]["molid"][molid]["coords"][:,2] -= traj_info[ts]["zlo"]
			bond_map = copy.deepcopy(traj_info[ts]["molid"][molid]["bond_map"])
			unwrap_molecule(traj_info[ts]["molid"][molid]["coords"], bond_map, box_dims, sr_llim)
			sr_llim += len(traj_info[ts]["molid"][molid]["coords"])
		print(f"Unwrapped!", flush=True)
	return

def coarse_grain_traj(traj_info):

	coarse_traj_info = copy.deepcopy(traj_info)
	for ts in traj_info:
		print (f"Coarse-graining timestep = {ts}...", end=' ', flush=True)
		dL = np.array([traj_info[ts]["xhi"]-traj_info[ts]["xlo"], traj_info[ts]["yhi"]-traj_info[ts]["ylo"], traj_info[ts]["zhi"]-traj_info[ts]["zlo"]], dtype=np.float64)
		for molid in traj_info[ts]["molid"]:
			coarse_traj_info[ts]["molid"][molid]["coords"]   = np.empty((0,3))
			coarse_traj_info[ts]["molid"][molid]["masses"]   = np.empty(0)
			coarse_traj_info[ts]["molid"][molid]["bond_map"] = {}
			coarse_traj_info[ts]["molid"][molid]["mol_name"] = traj_info[ts]["molid"][molid]["mol_name"]
			if molid == 1:
				# make the 19 monomers
				for i in range(30):
					if i==0:
						coords  = np.sum(traj_info[ts]["molid"][molid]["coords"][0:21]*traj_info[ts]["molid"][molid]["masses"][0:21][:,np.newaxis], axis=0)/np.sum(traj_info[ts]["molid"][molid]["masses"][0:21])
						coords %= dL
						coords[0] += traj_info[ts]["xlo"]
						coords[1] += traj_info[ts]["ylo"]
						coords[2] += traj_info[ts]["zlo"]
						coarse_traj_info[ts]["molid"][molid]["coords"]      = np.vstack((coarse_traj_info[ts]["molid"][molid]["coords"],coords))
						coarse_traj_info[ts]["molid"][molid]["masses"]      = np.hstack((coarse_traj_info[ts]["molid"][molid]["masses"], np.sum(traj_info[ts]["molid"][molid]["masses"][0:21])))
						coarse_traj_info[ts]["molid"][molid]["bond_map"][i+1] = [i+1+1]
					elif i<29:
						coords  = np.sum(traj_info[ts]["molid"][molid]["coords"][20+(i-1)*19:20+i*19]*traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:20+i*19][:,np.newaxis], axis=0)/np.sum(traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:20+i*19])
						coords %= dL
						coords[0] += traj_info[ts]["xlo"]
						coords[1] += traj_info[ts]["ylo"]
						coords[2] += traj_info[ts]["zlo"]
						coarse_traj_info[ts]["molid"][molid]["coords"]      = np.vstack((coarse_traj_info[ts]["molid"][molid]["coords"],coords))
						coarse_traj_info[ts]["molid"][molid]["masses"]      = np.hstack((coarse_traj_info[ts]["molid"][molid]["masses"], np.sum(traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:20+i*19])))
						coarse_traj_info[ts]["molid"][molid]["bond_map"][i+1] = [i+1-1,i+1+1]
					else:
						coords  = np.sum(traj_info[ts]["molid"][molid]["coords"][20+(i-1)*19:]*traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:][:,np.newaxis], axis=0)/np.sum(traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:])
						coords %= dL
						coords[0] += traj_info[ts]["xlo"]
						coords[1] += traj_info[ts]["ylo"]
						coords[2] += traj_info[ts]["zlo"]
						coarse_traj_info[ts]["molid"][molid]["coords"]      = np.vstack((coarse_traj_info[ts]["molid"][molid]["coords"],coords))
						coarse_traj_info[ts]["molid"][molid]["masses"]      = np.hstack((coarse_traj_info[ts]["molid"][molid]["masses"], np.sum(traj_info[ts]["molid"][molid]["masses"][20+(i-1)*19:])))
						coarse_traj_info[ts]["molid"][molid]["bond_map"][i+1] = [i+1-1]

				pass
				
			else:
				rcom = np.sum((traj_info[ts]["molid"][molid]["coords"]*traj_info[ts]["molid"][molid]["masses"][:,np.newaxis]), axis=0)/np.sum(traj_info[ts]["molid"][molid]["masses"])
				rcom %= dL
				rcom[0] += traj_info[ts]["xlo"]
				rcom[1] += traj_info[ts]["ylo"]
				rcom[2] += traj_info[ts]["zlo"]
				coarse_traj_info[ts]["molid"][molid]["coords"]   = rcom.reshape((1,-1))
				coarse_traj_info[ts]["molid"][molid]["bond_map"] = {}
				coarse_traj_info[ts]["molid"][molid]["masses"]   = np.sum(traj_info[ts]["molid"][molid]["masses"])
		print("Coarse grained!", flush=True)
	return coarse_traj_info

def write_cg_traj(coarse_traj_info, sim_data, cg_traj_filename, natoms):

	f = open(cg_traj_filename, 'w')

	for idx,ts in enumerate(coarse_traj_info):
		f.write("ITEM: TIMESTEP\n")
		f.write(f"{int(ts)}\n")
		f.write("ITEM: NUMBER OF ATOMS\n")
		f.write(f"{natoms}\n")
		f.write("ITEM: BOX BOUNDS pp pp pp\n")
		f.write(f'{coarse_traj_info[ts]["xlo"]}    {coarse_traj_info[ts]["xhi"]}\n')
		f.write(f'{coarse_traj_info[ts]["ylo"]}    {coarse_traj_info[ts]["yhi"]}\n')
		f.write(f'{coarse_traj_info[ts]["zlo"]}    {coarse_traj_info[ts]["zhi"]}\n')
		f.write("ITEM: ATOMS id type x y z\n")
		srno = 0
		for molid in coarse_traj_info[ts]["molid"]:
			for sr_atom, atom in enumerate(coarse_traj_info[ts]["molid"][molid]["coords"]):
				srno += 1
				f.write(f'{srno} {coarse_traj_info[ts]["molid"][molid]["mol_name"]} {atom[0]} {atom[1]} {atom[2]}\n')
	f.close()

	return

if __name__=="__main__":

	abs_start = time.time()
	start = time.time()
	print("Creating the topology object...", flush=True)
	if args.top == None:
		sim_info  = create_simulation_object(args.data)
		with open("top.pkl", 'wb') as f:
			pickle.dump(sim_info, f)
	else:
		with open(args.top, 'rb') as f:
			sim_info = pickle.load(f)

	print(sim_info)
	print(f"Time to make the topology is {time.time()-start} seconds.", flush=True)

	start = time.time()
	timestep_arr = traj_poker(args.traj)
	print(timestep_arr)
	print(f"Time taken to poke traj file is {time.time()-start} seconds.", flush=True)

	start = time.time()
	print("Creating the trajectory object...", flush=True)
	if args.trajpkl == None:
		traj_info  = create_traj_object(sim_info, args.traj, timestep_arr, args.nproc)
		with open("traj.pkl", 'wb') as f:
			pickle.dump(traj_info, f)
	else:
		with open(args.top, 'rb') as f:
			traj_info = pickle.load(f)
	print (f"traj_info = {traj_info}")
	print (f"Processed!.\nTime to make the traj object is {time.time()-start} seconds.", flush=True)

	print ("Created both the simulation object and traj object!", flush=True)

	start = time.time()
	unwrap_trajectory(traj_info)
	print(f"Time to unwrap trajectory is {time.time()-start} seconds.", flush=True)

	start = time.time()
	print(f"Writing the cg'd traj...", flush=True)
	cg_traj = coarse_grain_traj(traj_info)
	print(f"Coarse-grained the trajectory in {time.time()-start} seconds.", flush=True)

	start = time.time()
	print(f"Writing out cg'd traj...", flush=True)
	write_cg_traj(cg_traj, sim_info, args.ntraj, sim_info["nmolecules"]+29)
	print(f"Wrote out coarse-grained trajectory in {time.time()-start} seconds.", flush=True)

	print(f"Time to run computation was {time.time()-abs_start} seconds.", flush=True)

