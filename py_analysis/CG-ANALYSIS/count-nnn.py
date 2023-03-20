#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse
import copy
import sys
sys.path.insert (0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux
import os
import re
from numba import jit


parser = argparse.ArgumentParser(description="Run a coarse-grained simulation of the next-nearest neighbor variety.")
# parser.add_argument ("--tenergy-file", dest="energy", action='store', type=str, help="Provide address for energy dump file.")
parser.add_argument ("--tcoords-file", dest="coords", action='store', type=str, help="Provide address for coordinates file.")
# parser.add_argument ("--dump-loc", dest="dump", action="store", type=str, help="Where to dump next nearest coordinates.")
# parser.add_argument ("--model", dest='m', action='store', type=int, help="Select model directory.")
parser.add_argument ("--dop", dest='dop', action='store', type=int, help="Select degree of polymerization.")
# parser.add_argument ("-T", dest='T', action='store', type=float, help="Select temperature.")
# parser.add_argument ("-s", dest='s', action='store', type=int, help="Number of values to consider.")
args = parser.parse_args()


next_nearest_neighbor = np.array ([[2, 0, 0], [2, 1, 0], [2, -1, 0], [2, 0, 1], [2, 0, -1], [2, 1, 1], [2, -1, 1], [2, 1, -1], [2, -1, -1], \
[-2, 0, 0], [-2, 1, 0], [-2, -1, 0], [-2, 0, 1], [-2, 0, -1], [-2, 1, 1], [-2, -1, 1], [-2, 1, -1], [-2, -1, -1], [0, 2, 0], \
[1, 2, 0], [-1, 2, 0], [0, 2, 1], [0, 2, -1], [1, 2, 1], [-1, 2, 1], [1, 2, -1], [-1, 2, -1], [0, -2, 0], [1, -2, 0], [-1, -2, 0], \
[0, -2, 1], [0, -2, -1], [1, -2, 1], [-1, -2, 1], [1, -2, -1], [-1, -2, -1], [0, 0, 2], [0, 1, 2], [0, -1, 2], [1, 0, 2], [-1, 0, 2], \
[1, 1, 2], [-1, 1, 2], [1, -1, 2], [-1, -1, 2], [0, 0, -2], [0, 1, -2], [0, -1, -2], [1, 0, -2], [-1, 0, -2], [1, 1, -2], [-1, 1, -2], \
[1, -1, -2], [-1, -1, -2], [2, 2, 0], [2, 2, 1], [2, 2, -1], [-2, 2, 0], [-2, 2, 1], [-2, 2, -1], [2, -2, 0], [2, -2, 1], [2, -2, -1], \
[-2, -2, 0], [-2, -2, 1], [-2, -2, -1], [0, 2, 2], [1, 2, 2], [-1, 2, 2], [0, -2, 2], [1, -2, 2], [-1, -2, 2], [0, 2, -2], [1, 2, -2], \
[-1, 2, -2], [0, -2, -2], [1, -2, -2], [-1, -2, -2], [2, 0, 2], [2, 1, 2], [2, -1, 2], [-2, 0, 2], [-2, 1, 2], [-2, -1, 2], [2, 0, -2], \
[2, 1, -2], [2, -1, -2], [-2, 0, -2], [-2, 1, -2], [-2, -1, -2], [2, 2, 2], [-2, 2, 2], [2, -2, 2], [2, 2, -2], [-2, -2, 2], [-2, 2, -2], \
[2, -2, -2], [-2, -2, -2]] )

def get_starting_ind (filename, s):

	print ("energy file name = "+filename)
	df = pd.read_csv (filename, sep=' \| ', names=["energy", "mm_1", "mm_2", "ms1_tot", "time_step"], engine='python', skiprows=0)
	return int(df["time_step"].values[-s])


def infiltrate_coords_get_next_nearest_contacts (coords_file, starting_index, dop):

	edge     = aux.edge_length (dop)
	master_dict = aux.get_pdict (coords_file, starting_index, dop, edge, edge, edge)
	next_nearest_contacts = []

	for key in master_dict:
		polymer = master_dict[key][0]
		polymer = aux.unfuck_polymer (polymer, edge, edge, edge)
		# print ("key = ",key)
		# count = next_neighbor (polymer)
		count = 0
		for monomer in polymer:
		 	# print (monomer)
			for n_neigh in next_nearest_neighbor:
				next_neigh = monomer + n_neigh
				finder = len ( np.flatnonzero ( (polymer==next_neigh).all(1) ) )
				if finder == 1:
					count += 0.5
				elif finder > 1:
					print ("There is something fucked here.")
					exit ()
		next_nearest_contacts.append (count)

	return np.array (next_nearest_contacts)


if __name__=="__main__":

	dop = args.dop
	target_coords = args.coords
	
	# target next nearest neighbor contacts 
	next_neighbor_contacts_target = infiltrate_coords_get_next_nearest_contacts (target_coords, 0, dop)
	cdict = dict() 
	cdict ["timestep"] = np.arange(len(next_neighbor_contacts_target))
	cdict ["nnn"] = next_neighbor_contacts_target 
	cdf = pd.DataFrame.from_dict (cdict)
	print (cdf)
	cdf.to_csv ("nextn_dump_1.mc", index=False, header=False, sep='|')
	
