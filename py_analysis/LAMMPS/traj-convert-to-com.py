import numpy as np 
import re
import argparse

parser = argparse.ArgumentParser (description="Edit the data file to put the polymer in the center of the box.")
parser.add_argument ("--trajfile", dest='traj', action='store', type=str, help="Name of trajfile.")
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
				sim_data["molecules"][int(content[1])] = []
			sim_data["molecules"][int(content[1])].append(int(content[0]))

			continue
		elif extract_bonds:
			content = line.strip().split()
			sim_data["bonds"].append((int(content[2]), int(content[3])))
			continue

	sim_data["nmolecules"] = len(sim_data["molecules"])

	return sim_data

def create_traj_object(trajfile):

	traj = {}

	timestep_flag  = False
	numatoms_flag  = False
	boxbounds_flag = False
	itematoms_flag = False

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
			traj[ts]["coords"] = {}
			timestep_flag  = False
			numatoms_flag  = False
			boxbounds_flag = False
			itematoms_flag = True
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
			contents = line.strip().split()
			ts = int(contents[0])
			traj[ts] = {}
			continue

		elif numatoms_flag:
			contents = line.strip().split()
			traj[ts]["natoms"] = contents[0]
			continue

		elif itematoms_flag:
			contents = line.strip().split()
			traj[ts]["coords"][int(contents[0])] = {}
			traj[ts]["coords"][int(contents[0])]["id"]  = contents[1]
			traj[ts]["coords"][int(contents[0])]["loc"] = np.array([float(contents[2]), float(contents[3]), float(contents[4])], dtype=np.float64)
			continue

	return traj

def unwrap_molecule(molecule, bonds, box_dims):

	untested = list(range(len(polymer)))
	tested   = []
	queue    = []

	while untested:
		
		wait = []
		if not queue:
			queue.append(untested[0])

		for i in queue:
			neighbors = bonds[i]
			neighbors = [ni for ni in neighbors if ni not in tested]
			ri = polymer[i]
			for j in neighbors:
				rj = polymer[j]
				dr = rj - ri
				shift = np.round(dr/box_dims)
				polymer[j] -= shift*box_dims
			tested.append(i)
			untested.remove(i)
			wait.extend(neighbors)
		queue = list(set(wait[:]))

	return polymer 

if __name__=="__main__":

	print("Creating the simulation object...")
	sim_info  = create_simulation_object(args.data)
	print(sim_info)

	print("Creating the trajectory object...")
	traj_info = create_traj_object(args.traj)

	print ("Created both!")


