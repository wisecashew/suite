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
sys.path.insert (0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/py_analysis')
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


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument ("--model", dest='m', type=int, help="Select model directory.")
args = parser.parse_args()

if __name__=="__main__":

	info = aux.get_info ("TARGET/geom_and_esurf.txt")
	x = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]
	
	beta = 1.0
	chi = 0.05
	# obtain target parameters 
	energy_target = np.array  ( aux.get_energy ("TARGET/geom_and_esurf.txt") )
	energy_model  = np.array  ( aux.get_energy_cg2 ("MODEL"+str(args.m)+"/geom_and_esurf.txt") )
	energy_upd    = copy.copy ( energy_model )
	# energy_upd[2] = 0
	# get the target contacts 
	df_target = pd.read_csv ("TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "time_step"], skiprows=0)
	df_model  = pd.read_csv ("MODEL"+str(args.m)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "time_step"], skiprows=0)
	print ("Energetic parameters initial =", energy_upd)

	# now perform the update
	energy_upd[0] = energy_upd[0] - chi * ( np.mean ( df_target["mm_aligned"].values[-2000:] )  - np.mean ( df_model["mm_aligned"].values[-2000:] ) ) / ( beta * np.mean (df_model["mm_aligned"].values[-2000:]**2) - beta * np.mean(df_model["mm_aligned"].values[-2000:])**2  )


	energy_upd[1] = energy_upd[1] - chi * ( np.mean ( df_target["mm_naligned"].values[-2000:] ) - np.mean ( df_model["mm_naligned"].values[-2000:] ) ) / ( beta * np.mean (df_model["mm_naligned"].values[-2000:]**2) - beta * np.mean (df_model["mm_naligned"].values[-2000:])**2 )

	# print ("nalign num   =", ( np.mean ( df_target["mm_naligned"].values[-2000:] ) - np.mean ( df_model["mm_naligned"].values[-2000:] ) ) )

	# print ("nalign denom =",( beta * np.mean (df_model["mm_naligned"].values[-2000:]**2) - beta * np.mean (df_model["mm_naligned"].values[-2000:])**2  ))


	print ("Energetic parameters final =", energy_upd)
	file_new = open ("MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
	file_new.write (f"x = {x}\n")
	file_new.write (f"y = {y}\n")
	file_new.write (f"z = {z}\n")
	file_new.write (f"kT = {T}\n")
	file_new.write (f"Emm_a = {energy_upd[0]}\n")
	file_new.write (f"Emm_n = {energy_upd[1]}\n")
	file_new.write (f"Ems   = {energy_upd[2]}\n")
	file_new.write ("END OF FILE")
	file_new.close()


