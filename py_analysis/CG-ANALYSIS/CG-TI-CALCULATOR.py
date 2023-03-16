#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse
import sys
import copy 
sys.path.insert (0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/py_analysis')
import aux
import os
import re

parser = argparse.ArgumentParser(description="Set up a thermodynamic integration schema.")
parser.add_argument ('--model', dest='m', type=int, action='store', help='Enter a model iteration.')
args = parser.parse_args()

if __name__=="__main__":

	# obtain simulation parameters for each TI point 
	nlambda = 30
	W = np.polynomial.legendre.leggauss (nlambda)
	lambdas = W[0]
	weights = W[1]
	du_dl   = np.zeros (nlambda)

	print ("lambdas = ",lambdas)
	info = aux.get_info ("TARGET/geom_and_esurf.txt")
	x = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]
	beta = 1

	# get parameters from TARGET (HR)
	energy_target_list     = np.array(aux.get_energy ("TARGET/geom_and_esurf.txt"))
	print ("high res hamiltonian = ",energy_target_list)

	# get parameters from MODEL (CG)
	energy_model_list_raw  = aux.get_energy_cg1 ("MODEL" + str(args.m) +"/geom_and_esurf.txt")

	# print ("target energies = ",energy_target_list)
	energy_model_list = copy.copy(energy_target_list)
	energy_model_list[0] = energy_model_list_raw[0]; energy_model_list[1] = energy_model_list_raw[0];
	energy_model_list[2] = energy_model_list_raw[1]; energy_model_list[3] = energy_model_list_raw[1];
	energy_model_list = np.array (energy_model_list)
	print ("cg hamiltonian = ",energy_model_list)

	deltaU = (energy_model_list - energy_target_list)/2
	print ("deltaU = ",deltaU)

	for i in range (1, nlambda+1, 1):
		df_model = pd.read_csv ("ThermodynamicIntegration/MODEL"+str(args.m)+"_TI/ITER"+str(i)+"/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot", "ms1s2_aligned", "ms1s2_naligned", "time_step"], skiprows=0)
		du_dl[i-1] = np.mean ( deltaU[0]*df_model["mm_tot"].values[-2000:] + deltaU[2]*df_model["ms1_aligned"].values[-2000:] + deltaU[3]*df_model["ms1_naligned"].values[-2000:] )

	dF = np.sum (weights * du_dl)

	df_target = pd.read_csv ("TARGET/energydump_1.mc", sep='\|', engine='python', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot", "ms1s2_aligned", "ms1s2_naligned", "time_step"], skiprows=0)
	dU = np.mean (df_target["mm_aligned"].values[-2000:]*energy_model_list[0] + df_target["mm_naligned"].values[-2000:]*energy_model_list[1] + df_target["ms1_aligned"].values[-2000:]*energy_model_list[2] + df_target["ms1_naligned"].values[-2000:]*energy_model_list[3]) - np.mean (df_target["energy"].values[-2000:])

	Srel = (beta*dU - beta*dF)/32
	print  ("deltaU           = {}".format (dU)  )
	print  ("deltaF           = {}".format (dF)  )
	print  ("Relative entropy = {}".format (Srel))

	# plt.plot (lambdas, du_dl, color='steelblue', marker='o')
	# plt.ylabel ("$\\langle dU/d\\lambda \\rangle_{\\lambda}$")
	# plt.xlabel ("$\\lambda$")
	# plt.show()
