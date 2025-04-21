#!/home/satyend/.conda/envs/phase/bin/python

import lmp_coords
import numpy as np
import pandas as pd
import re
import os
import multiprocessing
import itertools
import time

import argparse 
parser = argparse.ArgumentParser(description="Take tmp coordinate files and make them into clean files.")
parser.add_argument ("--natom",    dest='n',   type=int, action='store', help="Number of atoms in the polymer.")
parser.add_argument ("--delta_ts", dest='d',   type=int, action='store', help="Enter difference between timesteps.")
parser.add_argument ("--dir",      dest='dir', type=str, action='store', help="Enter directory to dump cleaned coordinates.")
args = parser.parse_args ()

def clean_coords(coords_file, DELTA_TS, natoms):
	configs = lmp_coords.Configurations(coords_file, DELTA_TS)
	idx = coords_file.strip().split(".")[1]
	write_to_file = open(args.dir+"/polymer."+str(idx)+".lammpstrj", 'w')
	# print(f"Number of timesteps = {len(configs.timesteps)}", flush=True)
	for i in range(len(configs.timesteps)):
		write_to_file.write(f"{configs.r_timestep}\n")
		write_to_file.write(f"{configs.timesteps[i]}\n")
		write_to_file.write(f"{configs.r_natoms}\n")
		write_to_file.write(f"{natoms}\n")
		write_to_file.write(f"{configs.r_boxdims}\n")
		write_to_file.write(f"{configs.xdim[i][0]} {configs.xdim[i][1]}\n")
		write_to_file.write(f"{configs.ydim[i][0]} {configs.ydim[i][1]}\n")
		write_to_file.write(f"{configs.zdim[i][0]} {configs.zdim[i][1]}\n")
		write_to_file.write(f"{configs.r_colnames}\n")
		for j in range(len(configs.coords[i])):
			if j < natoms:
				write_to_file.write(f"{configs.coords[i][j][0]} {configs.coords[i][j][1]} {configs.coords[i][j][2]} {configs.coords[i][j][3]} {configs.coords[i][j][4]}\n")
			else:
				break
	write_to_file.close()
	return

if __name__=="__main__":
	
	start = time.time()

	# get the list of files
	files = os.listdir(".")

	# get all the trajectories
	coords_files = []
	for file in files:
		if re.match(r"^tmp\.[0-9]+\.lammpstrj", file):
			coords_files.append(file)

	# sort the trajectories
	coords_files.sort(key=lambda x:int(x.split(".")[1]))

	# print out the coordinate files
	print(f"Number of coords files to process is {len(coords_files)}.", flush=True)

	# get the frequency of dumps
	DELTA_TS = args.d

	# define the number of atoms in the polymer
	natoms = args.n

	# parse through the files
	# for idx, coords_file in enumerate(coords_files[0:1]):
	for idx in range(len(coords_files)):
		print(f"At index {idx}.", flush=True)
		clean_coords(coords_files[idx], DELTA_TS, natoms)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

