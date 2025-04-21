#!/home/satyend/.conda/envs/phase/bin/python

import lmp_data
import lmp_settings
import MDAnalysis as mda
import freesasa
import numpy as np
import argparse
import time

parser = argparse.ArgumentParser(description="Compute the solvent accessible surface area.")
parser.add_argument("--datafile",     dest="df",    type=str, action='store', help="enter address of datafile.")
parser.add_argument("--settingsfile", dest="sf",    type=str, action='store', help="enter address of settingsfile.")
parser.add_argument("--coords",       dest='c',     type=str, action='store', help="enter address of coordinate file.")
parser.add_argument("--sasadump",     dest="sd",    type=str, action='store', help="enter address of SASA dump file.")
args = parser.parse_args()

if __name__=="__main__":
	
	start = time.time()

	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	my_settings = lmp_settings.Settings(settings)
	my_settings.read_settings()

	my_data     = lmp_data.Data(topo)
	my_data.read_data()

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	print("Begin selecting some atoms for the polymer...", end=' ', flush=True)
	polymer = u.select_atoms('resid 1')
	print("done!")

	# positions in MDA are shifted so that xlo, ylo, zlo = 0
	# shift = np.array([1.5620737849713620e-01, 1.5620737849713620e-01, 1.5620737849713620e-01])
	vdw_radii = list()
	for atom in polymer.atoms:
		vdw_radii.append(my_settings.pair_coeff[int(atom.type)][1])

	parameters = freesasa.Parameters()

	print(f"Algorithm is {parameters.algorithm()}.")
	print(f"Probe radius is {parameters.probeRadius()}Å.")
	
	sasa_container = list()
	timesteps      = list()

	for idx, ts in enumerate(u.trajectory):
		print(f"timestep = {ts.time/2} | time = {ts.time} nanoseconds.", flush=True)
		polymer.unwrap()
		coords = polymer.positions 
		coords = coords.flatten().tolist()
		result = freesasa.calcCoord(coords, vdw_radii)
		timesteps.append(ts.time/2)
		sasa_container.append(result.totalArea())

	dumpfile = open(args.sd, 'w')
	for i in range(len(sasa_container)):
		if i < len(sasa_container)-1:
			dumpfile.write("{:.2f} {:.2f}\n".format(timesteps[i], sasa_container[i]))
		else:
			dumpfile.write("{:.2f} {:.2f}".format(timesteps[i], sasa_container[i]))
	dumpfile.close()
	# print(sasa_container)
	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)
