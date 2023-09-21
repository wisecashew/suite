#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-N", dest='N', type=int, action='store', help="Length of polymer.")
parser.add_argument ("--nskip", dest='nskip', type=int, action='store', help="Number of timesteps to skip.")
args = parser.parse_args ()

def f_exact (n_, t_, k_, b_, p_, x_):

	t1 = 2 * np.log (x_) 
	t2 = - 6 * (np.sin(p_*np.pi/(2*n_))**2) / (b_ ** 2) * (x_ ** 2)
	t3 = np.log(24 * np.sqrt(6/np.pi) * (np.sin (np.pi*p_/(2*n_))/b_)  ** 3)

	return -k_*t_*(t1+t2+t3)


def Xp (traj, p, N):
	xp = []
	i_list = np.arange (0, N)
	print (i_list + 1/2)
	cos_postfix = np.cos (p*np.pi/N * (i_list+1/2))
	print (f"N = {N}")
	prefactor = np.sqrt(2/N)

	keys_to_probe = list (traj.keys())[1000:]

	for key in keys_to_probe:
		if len(traj[key])!=21:
			print (f"traj[key].shape = {traj[key].shape}")
			print ("Something's fucked.")
			exit ()

		x1_vec = np.sum (prefactor * traj[key] * cos_postfix[:, np.newaxis], axis=0)
		x1_v_mag = np.linalg.norm (x1_vec)
		xp.append ( np.linalg.norm (x1_v_mag) )

	return xp

def average_successive_row_distance(array):
    # Check if the input array is empty or has only one row
    if len(array) < 2:
        raise ValueError("Input array must have at least two rows.")

    # Calculate the pairwise Euclidean distances between successive rows
    row_diff = np.diff(array, axis=0)
    distances = np.square(np.linalg.norm(row_diff, axis=1))

    # Compute the average distance
    average_distance = np.mean(distances)

    return average_distance

def bond_l (traj, nskip):
	timesteps = list(traj.keys())
	bondls = []
	for ts in timesteps[nskip:]:
		print (f"@ timestep = {ts}...")
		bondls.append(average_successive_row_distance(traj[ts]))

	avg_bonds = np.mean(bondls)
	return avg_bonds

if __name__=="__main__":

	f = open ("coords.nvt.lammpstrj", 'r')
	coords_flag    = False
	timestep_flag  = False
	trajectory = {}

	DOP = args.N

	for line in f:
		
		if re.findall ("^ITEM: TIMESTEP$", line):
			polymer = np.zeros ((DOP,3))
			timestep_flag     = True
			coords_flag       = False
			continue
		
		elif re.findall ("^[0-9]+$", line) and timestep_flag:
			timestep = int (line)
			continue
		
		elif re.findall ("^ITEM: NUMBER OF ATOMS$", line):
			timestep_flag = False
			continue

		elif re.findall ("^ITEM: ATOMS id type x y z$", line):
			coords_flag   = True
			continue
		
		elif re.findall ("^.*\\s.*\\s.*\\s.*\\s.*$", line) and coords_flag:
			info = line.strip().split()
			polymer[int(info[0])-1, :] = np.array ([float(info[2]), float(info[3]), float(info[4])])
			if int(info[0]) == DOP:
				trajectory [timestep] = polymer
			continue
		
		else:
			continue

	avg_bondl = bond_l(trajectory, args.nskip)
	print (f"average bond length = {np.sqrt(avg_bondl)}.")
