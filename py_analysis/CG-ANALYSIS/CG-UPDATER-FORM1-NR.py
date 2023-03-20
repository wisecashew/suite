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

parser = argparse.ArgumentParser(description="Run a coarse-grained simulation with isotropic simulation.")
parser.add_argument ("--model", dest='m', type=int, help="Select model directory.")
parser.add_argument ("-T", dest='T', type=float, help="Select temperature.")
args = parser.parse_args()

if __name__=="__main__":

	info = aux.get_info (str(T)+"/TARGET/geom_and_esurf.txt")
	x = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]
	k    = 1
	T    = args.T
	beta = 1/(k*T)
	chi  = 0.05

	# obtain target parameters
	energy_target = np.array  ( aux.get_energy_target (str(T)+"/TARGET/geom_and_esurf.txt") )

	# obtain model parameters
	energy_model  = np.array  ( aux.get_energy_form_2 (str(T)+"/FORM1/MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )

	# get the target contacts 
	df_target = pd.read_csv (str(T)+"/TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned","ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], skiprows=0)
	df_model  = pd.read_csv (str(T)+"/FORM1/MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "ms1_tot", "time_step"], skiprows=0)

	print ("Energetic parameters initial =", energy_upd)
	diff_mm = np.mean (df_model["mm_tot"].values[-2000:]) - np.mean (df_target["mm_tot"].values[-2000:])
	delta_file = open (str(T)+"/FORM1/MODEL"+str(args.m+1)+"/delta.mc", 'w')
	delta_file.write  (f"<N_mm>_model - <N_mm>_target = {diff_mm}")
	delta_file.close  ()

	# now perform the update
	energy_upd[0] = energy_upd[0] - chi * ( np.mean ( df_target["mm_tot"].values[-2000:] )  - np.mean ( df_model["mm_tot"].values[-2000:] ) ) / ( beta * np.mean (df_model["mm_tot"].values[-2000:]**2) - beta * np.mean(df_model["mm_tot"].values[-2000:])**2  )


	print ("Energetic parameters final =", energy_upd)
	file_new = open (str(T)+"/FORM1/MODEL"+str(args.m+1)+"/geom_and_esurf_new.txt", 'w')
	file_new.write (f"x = {x}\n")
	file_new.write (f"y = {y}\n")
	file_new.write (f"z = {z}\n")
	file_new.write (f"kT = {T}\n")
	file_new.write (f"Emm   = {energy_upd[0]}\n")
	file_new.write (f"Ems   = 0\n")
	file_new.write ("END OF FILE")
	file_new.close()


