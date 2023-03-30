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

def delta_scale_finder (filename):
	with open (filename) as f:
		for line in f:
			if re.search ("<N_mm>", line):
				delta = re.search ("-\d+\.\d+|\d+\.\d+", line)
				delta = float(delta.group(0))
			elif re.search ("scale", line):
				scale = re.search ("\d+\.\d+|\d+", line)
				scale = float (scale.group(0))
	return delta, scale


parser = argparse.ArgumentParser(description="Run a coarse-grained simulation with directional interactions.")
parser.add_argument ("--model", dest='m', action='store', type=int, help="Select model directory.")
parser.add_argument ("-T", dest='T', action='store', type=float, help="Select temperature.")
parser.add_argument ("-s", dest='s', action='store', type=int, help="Number of values to consider.")
args = parser.parse_args()

def det_chi (T, error_scale):

	if error_scale == 0.0:
		return 0
	print ("temperature = ",T)
	print ("es = ", error_scale)
	log_T     = np.ceil(np.log10(T))
	log_error = np.ceil (np.log10(error_scale))

	if log_T >= 0:
		chi      = 5 * (10 ** (log_beta - 2) )
	else:
		chi = 5 * (10 ** (log_beta - 1))

	if log_error > log_beta+1:
		chi *= 1.1
	elif log_error < log_beta:
		chi *= 0.75
	
	return chi

if __name__=="__main__":

	k           = 1
	T           = args.T
	s           = args.s
	info        = aux.get_info (str(T)+"/TARGET/geom_and_esurf.txt")
	x           = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]
	beta        = 1/(k*T)

	# obtain target parameters 
	energy_target = np.array  ( aux.get_energy_target (str(T)+"/TARGET/geom_and_esurf.txt") )
	energy_model  = np.array  ( aux.get_energy_form_2 (str(T)+"/FORM2/MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )

	# get the target contacts 
	num_list_T = aux.dir2nsim (os.listdir(str(T)+"/TARGET/."))
	N_mm_target_list  = []
	N_mma_target_list = []
	N_mmn_target_list = []
	for i in num_list_T:
		df_target = pd.read_csv (str(T)+"/TARGET/energydump_"+str(i)+".mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], skiprows=0)
		N_mm_target_list.extend ( list (df_target["mm_tot"].values[-s:]) )
		N_mma_target_list.extend( list (df_target["mm_aligned"].values[-s:]) )
		N_mmn_target_list.extend( list (df_target["mm_naligned"].values[-s:]) )

	N_mm_target_list  = np.array (N_mm_target_list)
	N_mma_target_list = np.array (N_mma_target_list)
	N_mmn_target_list = np.array (N_mmn_target_list)

	avg_N_mm_target  = np.mean(N_mm_target_list)
	avg_N_mma_target = np.mean(N_mma_target_list)
	avg_N_mmn_target = np.mean(N_mmn_target_list)

	num_list_M = aux.dir2nsim (os.listdir (str(T)+"/FORM2/MODEL"+str(args.m)+"/."))
	N_mm_model_list  = []
	N_mma_model_list = []
	N_mmn_model_list = []

	for i in num_list_M:
		df_model  = pd.read_csv ( str(T)+"/FORM2/MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "time_step"], skiprows=0)
		N_mm_model_list.extend  ( list(df_model["mm_tot"].values[-s:]) )
		N_mma_model_list.extend ( list(df_model["mm_aligned"].values[-s:]) )
		N_mmn_model_list.extend ( list(df_model["mm_naligned"].values[-s:]) )

	N_mm_model_list  = np.array(N_mm_model_list)
	N_mma_model_list = np.array(N_mma_model_list)
	N_mmn_model_list = np.array(N_mmn_model_list)

	avg_N_mm_model  = np.mean(N_mm_model_list)
	avg_N_mma_model = np.mean(N_mma_model_list)
	avg_N_mmn_model = np.mean(N_mmn_model_list)
	chi_aligned  = 0.1
	chi_naligned = 0.1

	diff_mm  = (avg_N_mm_model  - avg_N_mm_target  )
	diff_mma = (avg_N_mma_model - avg_N_mma_target )
	diff_mmn = (avg_N_mmn_model - avg_N_mmn_target )

	delta_file = open (str(T)+"/FORM2/MODEL"+str(args.m)+"/delta.mc", 'w')
	delta_file.write  (f"beta = {beta}\n")
	delta_file.write  (f"<N_mm>_model   - <N_mm>_target   = {diff_mm}\n" )
	delta_file.write  (f"<N_mm_a>_model - <N_mm_a>_target = {diff_mma}\n")
	delta_file.write  (f"<N_mm_n>_model - <N_mm_n>_target = {diff_mmn}\n")

	# now perform the update
	aligned_num    = avg_N_mma_target  - avg_N_mma_model
	aligned_denom  = beta * np.mean (N_mma_model_list ** 2) - beta * np.mean (N_mma_model_list) ** 2

	if np.abs(aligned_denom) < 1e-4:
		aligned_denom = 0.01
	print ("aligned num =", aligned_num)
	energy_upd[0] = energy_upd[0] - chi_aligned  * aligned_num / aligned_denom

	naligned_num   = avg_N_mmn_target - avg_N_mmn_model
	naligned_denom = beta * np.mean (N_mmn_model_list ** 2) - beta * np.mean (N_mmn_model_list) ** 2

	if np.abs (naligned_denom) < 1e-4:
		naligned_denom = 0.01

	energy_upd[1] = energy_upd[1] - chi_naligned * naligned_num / naligned_denom

	delta_file.write ("aligned update numerator   = {}\n".format (aligned_num) )
	delta_file.write ("aligned update denominator = {}\n".format (aligned_denom) )
	delta_file.write ("aligned chi    = {}\n".format (chi_aligned))
	delta_file.write ("aligned update = {}\n".format (chi_aligned  * aligned_num / aligned_denom) )

	delta_file.write ("naligned update numerator   = {}\n".format (naligned_num) )
	delta_file.write ("naligned update denominator = {}\n".format (naligned_denom) )
	delta_file.write ("naligned chi    = {}\n".format (chi_naligned) )
	delta_file.write ("naligned update = {}\n".format (chi_naligned  * naligned_num / naligned_denom) )
	delta_file.close ()

	print ("Energetic parameters final =", energy_upd)

	file_new = open (str(T)+"/FORM2/MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
	file_new.write (f"x = {x}\n" )
	file_new.write (f"y = {y}\n" )
	file_new.write (f"z = {z}\n" )
	file_new.write (f"kT = {T}\n")
	file_new.write ("Emm_a = {:2.10f}\n".format(energy_upd[0]))
	file_new.write ("Emm_n = {:2.10f}\n".format(energy_upd[1]))
	file_new.write (f"Ems   = 0\n")
	file_new.write ("END OF FILE" )
	file_new.close()


