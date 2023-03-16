#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
	print ("lambdas = ",lambdas)
	info = aux.get_info ("TARGET/geom_and_esurf.txt")
	x = info[0]; y = info[1]; z = info[2]; T = info[3]; frac = info[4]

	# get parameters from TARGET
	energy_target_list     = np.array(aux.get_energy ("TARGET/geom_and_esurf.txt"))
	energy_model_list_raw  = aux.get_energy_cg ("MODEL" + str(args.m) +"/geom_and_esurf.txt")
	print ("target energies = ",energy_target_list)
	energy_model_list = copy.copy(energy_target_list)
	energy_model_list[0] = energy_model_list_raw[0]; energy_model_list[1] = energy_model_list_raw[0];
	energy_model_list[2] = energy_model_list_raw[1]; energy_model_list[3] = energy_model_list_raw[1];
	print ("cg energies = ",energy_model_list)
	energy_model_list = np.array (energy_model_list)

	for i in range (nlambda):
		q_pts = (1-lambdas[i])/2 * energy_target_list + (1+lambdas[i])/2 * energy_model_list
		# print ("q_pts = ",q_pts)
		f = open ("ThermodynamicIntegration/MODEL" + str(args.m) + "_TI/ITER"+str(i+1)+"/geom_and_esurf.txt", 'w')
		f.write (f"x = {x}\n")
		f.write (f"y = {y}\n")
		f.write (f"z = {z}\n")
		f.write (f"frac = {frac}\n")
		f.write (f"kT = {T}\n")
		f.write ("m1-m1:isotropic\n")
		f.write ("m1-s1:parallel\n" )
		f.write ("m1-s2:isotropic\n")
		f.write ("s1-s2:isotropic\n")
		f.write (f"Emm_a = {q_pts[0]}\n")
		f.write (f"Emm_n = {q_pts[1]}\n")
		f.write (f"Ems1_a = {q_pts[2]}\n")
		f.write (f"Ems1_n = {q_pts[3]}\n")
		f.write (f"Ems2_a = {q_pts[4]}\n")
		f.write (f"Ems2_n = {q_pts[5]}\n")
		f.write (f"Es1s2_a = {q_pts[6]}\n")
		f.write (f"Es1s2_n = {q_pts[7]}\n")
		f.write (f"END OF FILE")
		f.close ()

