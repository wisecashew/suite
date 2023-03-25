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

	regularizer = 0.0
	k           = 1
	T           = args.T
	s           = args.s
	info        = aux.get_info (str(T)+"/TARGET/geom_and_esurf.txt")
	x           = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]
	beta        = 1/(k*T)

	# obtain target parameters 
	energy_target = np.array  ( aux.get_energy_target (str(T)+"/TARGET/geom_and_esurf.txt") )
	energy_model  = np.array  ( aux.get_energy_form_4 (str(T)+"/FORM4/MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )

	# energy_upd[2] = 0
	# get the target contacts 
	df_target = pd.read_csv (str(T)+"/TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], skiprows=0)
	df_target_magnet = pd.read_csv (str(T)+"/TARGET/magnetization_dump_1.mc", sep='\|', engine='python', names=["time", "magnetization"], skiprows=1)
	df_model  = pd.read_csv ( str(T)+"/FORM2/MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "magnetization", "time_step"], skiprows=0)

	aligned_error_scale  = np.abs(np.mean(df_target["mm_aligned" ].values[-s:]  - df_model["mm_aligned" ].values[-s:] ))
	naligned_error_scale = np.abs(np.mean(df_target["mm_naligned"].values[-s:]  - df_model["mm_naligned"].values[-s:] ))
	

	chi_aligned  = 0.1  # det_chi  (beta, aligned_error_scale)
	print (f"chi_aligned = {chi_aligned}")
	chi_naligned = 0.1 # det_chi (beta, naligned_error_scale)
	print (f"chi_naligned = {chi_naligned}")

	print ("Energetic parameters initial =", energy_upd)
	diff_mm  = -(np.mean (df_model["mm_tot"].values[-s:]) - np.mean (df_target["mm_tot"].values[-s:]))
	diff_mma = -(np.mean (df_model["mm_aligned"].values[-s:] ) - np.mean (df_target["mm_aligned"].values[-s:] ))
	diff_mmn = -(np.mean (df_model["mm_naligned"].values[-s:]) - np.mean (df_target["mm_naligned"].values[-s:]))
	diff_mag = -(np.mean (df_model["magnetization"].values[-s:]) - np.mean (df_target_magnet["magnetization"].values[-s:]))
	delta_file = open (str(T)+"/FORM2/MODEL"+str(args.m)+"/delta.mc", 'w')
	delta_file.write  (f"beta = {beta}\n")
	delta_file.write  (f"<N_mm>_model   - <N_mm>_target   = {diff_mm}\n" )
	delta_file.write  (f"<N_mm_a>_model - <N_mm_a>_target = {diff_mma}\n")
	delta_file.write  (f"<N_mm_n>_model - <N_mm_n>_target = {diff_mmn}\n")
	delta_file.write  (f"<M>_model      - <M>_target      = {diff_mag}\n")

	# now perform the update
	aligned_num    = np.mean ( df_target["mm_aligned"].values[-s:] )  - np.mean ( df_model["mm_aligned"].values[-s:] )
	aligned_denom  = beta * np.mean (df_model["mm_aligned"].values[-s:]**2) - beta * np.mean(df_model["mm_aligned"].values[-s:])**2

	if np.abs(aligned_denom) < 1e-4:
		aligned_denom = 0.01

	energy_upd[0] = energy_upd[0] - chi_aligned  * aligned_num / aligned_denom

	naligned_num   = np.mean ( df_target["mm_naligned"].values[-s:] ) - np.mean ( df_model["mm_naligned"].values[-s:] )
	naligned_denom = beta * np.mean (df_model["mm_naligned"].values[-s:]**2) - beta * np.mean (df_model["mm_naligned"].values[-s:])**2

	if np.abs (naligned_denom) < 1e-4:
		naligned_denom = 0.01

	energy_upd[1] = energy_upd[1] - chi_naligned * naligned_num / naligned_denom

	bfield_num    = np.mean (df_target_magnet["magnetization"].values[-s:]) - np.mean (df_model["magnetization"].values[-s:])
	bfield_den    = beta * np.mean (df_model["magnetization"].values[-s:]**2) - beta * np.mean (df_model["magnetization"].values[-s:])**2

	energy_upd[2] = energy_upd[2] - 0.1 * bfield_num / bfield_den

	delta_file.write ("aligned update numerator   = {}\n".format (aligned_num) )
	delta_file.write ("aligned update denominator = {}\n".format (aligned_denom ) )
	delta_file.write ("aligned chi    = {}\n".format (chi_aligned))
	delta_file.write ("aligned update = {}\n".format (-0.1 * aligned_num/aligned_denom) )

	delta_file.write ("naligned update numerator   = {}\n".format (naligned_num ) )
	delta_file.write ("naligned update denominator = {}\n".format (naligned_denom ) )
	delta_file.write ("naligned chi    = {}\n".format (chi_naligned) )
	delta_file.write ("naligned update = {}\n".format (-0.1 * naligned_num/naligned_denom ) )

	delta_file.write ("bfield update numerator = {}\n".format   (bfield_num))
	delta_file.write ("bfield update denominator = {}\n".format (bfield_den))
	delta_file.write ("bfield update = {}\n".format (-0.1 * bfield_num/bfield_den) )

	delta_file.close ()

	print ("Energetic parameters final =", energy_upd)

	file_new = open (str(T)+"/FORM2/MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
	file_new.write (f"x = {x}\n" )
	file_new.write (f"y = {y}\n" )
	file_new.write (f"z = {z}\n" )
	file_new.write (f"kT = {T}\n")
	file_new.write ("Emm_a = {:2.10f}\n".format(energy_upd[0]))
	file_new.write ("Emm_n = {:2.10f}\n".format(energy_upd[1]))
	file_new.write ("B     = {:2.10f}\n".format(energy_upd[2]))
	file_new.write ("END OF FILE" )
	file_new.close()


