#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import re
import argparse
import sys
sys.path.insert (0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/py_analysis')
import aux
import os


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument ("--model", dest='m', type=int, help="Select model directory.")
args = parser.parse_args()


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



if __name__=="__main__":

	# get the entire list of potential energy surfaces. 
	# step 1: find <dU_M/dlambdai>_T. 
	# T is the target ensemble: the real ensemble. 
	# M is the coarse-grained potential energy function. 
	# two lambdas: E_mm, E_ms. 
	# step 1 involves finding <N_mm> and <N_ms>
	# step 1 (recast): find <N_mm> and <N_ms> 
	temps = aux.dir2float ( os.listdir(".") )
	for T in temps:
		print (f"In T = {T}...", flush=True)
		df_T = pd.read_csv (str(T)+"/TARGET/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
		avg_Nmm_T = np.mean  (df_T["mm_tot" ].values[-2000:])
		print  ("target average mm contacts = " + str(avg_Nmm_T))
	
		print (f"\tIn MODEL{args.m}...", flush=True)
		df_M  = pd.read_csv (str(T)+"/MODEL"+str(args.m)+"/energydump_1.mc",  sep=' \| ', names=["energy", "mm_tot", "ms_tot", "time_step"], engine='python', skiprows=0)
		avg_Nmm_M = np.mean (df_M["mm_tot"].values[-2000:])

		delta = avg_Nmm_T - avg_Nmm_M
		f = open (str(T)+"/MODEL"+str(args.m)+"/delta.mc", 'w')
		f.write ("<N_mm>_T - <N_mm>_M = {}\n".format (delta))

		try:
			delta_old, scale_old = delta_scale_finder (str(T)+"/MODEL"+str(args.m-1)+"/delta.mc")
			if np.sign (delta) * np.sign (delta_old) == -1:
				scale = scale_old*0.5
			else:
				scale = scale_old 
		except FileNotFoundError:
			scale = 0.5
			pass

		f.write ("scale = {:2.10f}".format (scale))
		f.close()

