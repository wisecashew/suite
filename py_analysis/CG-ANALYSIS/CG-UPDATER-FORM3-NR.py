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


parser = argparse.ArgumentParser(description="Run a coarse-grained simulation of the next-nearest neighbor variety.")
parser.add_argument ("--model", dest='m', action='store', type=int, help="Select model directory.")
parser.add_argument ("-T", dest='T', action='store', type=float, help="Select temperature.")
parser.add_argument ("-s", dest='s', action='store', type=int, help="Number of values to consider.")
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

def get_starting_ind (temp, s):
    filename = str(temp) + "/TARGET/energydump_1.mc"
    df = pd.read_csv (filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    L = len(df["energy"])

    return int(df["time_step"].values[L-2000])


def infiltrate_coords_get_next_nearest_contacts (dop, temp, coords_file, starting_index):
	filename = str(temp) + "/TARGET/coords_1.mc"
	edge     = aux.edge_length (dop)
	master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)

	next_nearest_contacts = []

	for key in master_dict:
		polymer = master_dict[key][0]
		polymer = unfuck_polymer (polymer, edge, edge, edge)
		count = 0
		for monomer in polymer:
			for n_neigh in next_nearest_neighbor:
				next_neigh = monomer + n_neigh 
				if next_neigh  





if __name__=="__main__":

	T = args.T
	k = 1
	s = args.s
	info = aux.get_info (str(T)+"/TARGET/geom_and_esurf.txt")
	x = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]

	beta = 1/(k*T)
	chi = 0.05
	# obtain target parameters 
	energy_target = np.array  ( aux.get_energy (str(T)+"/TARGET/geom_and_esurf.txt") )
	energy_model  = np.array  ( aux.get_energy_form_3 (str(T)+"/FORM2/MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )

	# energy_upd[2] = 0
	# get the target contacts 
	df_target = pd.read_csv (str(T)+"/TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], skiprows=0)
	df_model  = pd.read_csv (str(T)"/FORM3/MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "time_step"], skiprows=0)


	print ("Energetic parameters initial =", energy_upd)
	diff_mm = np.mean (df_model["mm_tot"].values[-2000:]) - np.mean (df_target["mm_tot"].values[-2000:])

	delta_file = open (str(T)+"/FORM3/MODEL"+str(args.m+1)+"/delta_new.mc", 'w')
	delta_file.write (f"<N_mm>_model - <N_mm>_target = {diff_mm}")
	delta_file.close()

	# now perform the update
	energy_upd[0] = energy_upd[0] - chi * ( np.mean ( df_target["mm_aligned"].values[-2000:] )  - np.mean ( df_model["mm_aligned"].values[-2000:] ) ) / ( beta * np.mean (df_model["mm_aligned"].values[-2000:]**2) - beta * np.mean(df_model["mm_aligned"].values[-2000:])**2  )

	# read the target coords file and get the next nearest neighbor distribution
	

	# energy_upd[1] = energy_upd[1] - chi * ( np.mean ( df_target["mm_naligned"].values[-2000:] ) - np.mean ( df_model["mm_naligned"].values[-2000:] ) ) / ( beta * np.mean (df_model["mm_naligned"].values[-2000:]**2) - beta * np.mean (df_model["mm_naligned"].values[-2000:])**2 )

	# print ("nalign num   =", ( np.mean ( df_target["mm_naligned"].values[-2000:] ) - np.mean ( df_model["mm_naligned"].values[-2000:] ) ) )
	# print ("nalign denom =",( beta * np.mean (df_model["mm_naligned"].values[-2000:]**2) - beta * np.mean (df_model["mm_naligned"].values[-2000:])**2  ))

	print ("Energetic parameters final =", energy_upd)
	file_new = open (str(T)+"/FORM3/MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
	file_new.write (f"x = {x}\n")
	file_new.write (f"y = {y}\n")
	file_new.write (f"z = {z}\n")
	file_new.write (f"kT = {T}\n")
	file_new.write (f"Emm_1 = {energy_upd[0]}\n")
	file_new.write (f"Emm_2 = {energy_upd[1]}\n")
	file_new.write (f"Ems   = {energy_upd[2]}\n")
	file_new.write ("END OF FILE")
	file_new.close()


