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
# parser.add_argument ("--tenergy-file", dest="energy", action='store', type=str, help="Provide address for energy dump file.")
parser.add_argument ("--model", dest='m', action='store', type=int, help="Select model directory.")
parser.add_argument ("--dop", dest='dop', action='store', type=int, help="Select degree of polymerization.")
parser.add_argument ("-T", dest='T', action='store', type=float, help="Select temperature.")
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

	# filename = str(temp) + "/TARGET/energydump_1.mc"
	df = pd.read_csv (filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	return int(df["time_step"].values[-s])


def det_chi (beta, error_scale):
	log_beta = np.ceil(np.log10(beta))
	if log_beta >= 0:
		chi      = 5 * (10 ** (log_beta - 2))
	else:
		chi = 5 * (10 ** (log_beta - 1))
	log_error = np.ceil (np.log10(error_scale))
	if log_error > log_beta+1:
		chi *= 1.1
	elif log_error < log_beta:
		chi *= 0.75
	
	return chi



if __name__=="__main__":

	regularizer = 0.00
	dop         = 32
	T           = args.T
	k           = 1
	s           = 2000
	info        = aux.get_info (str(T)+"/TARGET/geom_and_esurf.txt")

	x    = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4];

	beta = 1/(k*T)
	chi = 0.1

	# obtain target parameters 
	energy_target = np.array  ( aux.get_energy_target (str(T)+"/TARGET/geom_and_esurf.txt") )
	energy_model  = np.array  ( aux.get_energy_form_3 (str(T)+"/FORM3/MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )

	# energy_upd[2] = 0
	# get the target contacts 
	df_target_n = pd.read_csv (str(T)+"/TARGET/nextn_dump_1.mc", sep='\|', engine='python', names=["time_step", "next_neighbor"], skiprows=0)
	df_target   = pd.read_csv (str(T)+"/TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], skiprows=0)
	df_model    = pd.read_csv (str(T)+"/FORM3/MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_1", "mm_2", "ms1_tot", "time_step"], skiprows=0)

	chi_nearest      = 0.1
	chi_next_nearest = 0.1

	print (f"chi_nearest = {chi_nearest}")
	print (f"chi_next_nearest = {chi_next_nearest}")

	neighbor_contacts_target      = df_target   ["mm_aligned"].values[-s:]
	neighbor_contacts_model       = df_model    ["mm_1"].values[-s:]
	next_neighbor_contacts_target = df_target_n ["next_neighbor"].values[-s:]
	next_neighbor_contacts_model  = df_model    ["mm_2"].values[-s:]

	print ("Energetic parameters initial =", energy_upd)
	diff_next = np.mean (df_model["mm_1"].values[-s:]) - np.mean (df_target["mm_tot"].values[-s:])
	diff_nearest_next  = -(np.mean(df_target_n["next_neighbor"].values[-s:]  - df_model["mm_2"].values[-s:]))

	delta_file = open (str(T)+"/FORM3/MODEL"+str(args.m)+"/delta.mc", 'w')
	delta_file.write (f"beta = {beta}\n")
	delta_file.write (f"<N_mm>_model - <N_mm>_target = {diff_next}\n")
	delta_file.write (f"<N_mm_n>_model - <N_mm_n>_target = {diff_nearest_next}\n")

	energy_upd[0] = energy_upd[0] - chi_nearest      * beta * (diff_next)
	energy_upd[1] = energy_upd[1] - chi_next_nearest * beta * (diff_nearest_next)

	delta_file.write (f"chi_nearest = {chi_nearest}\n")
	delta_file.write (f"nearest update = {energy_upd[0]}\n")
	delta_file.write (f"chi_next_nearest = {chi_next_nearest}\n")
	delta_file.write (f"next nearest update = {energy_upd[1]}\n")

	delta_file.close ()
	print ("Energetic parameters initial =", energy_upd)

	file_new = open (str(T)+"/FORM3/MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
	file_new.write (f"x  = {x}\n")
	file_new.write (f"y  = {y}\n")
	file_new.write (f"z  = {z}\n")
	file_new.write (f"kT = {T}\n")
	file_new.write ("Emm_1 = {:2.10f}\n".format(energy_upd[0]))
	file_new.write ("Emm_2 = {:2.10f}\n".format(energy_upd[1]))
	file_new.write ("END OF FILE")
	file_new.close ()


