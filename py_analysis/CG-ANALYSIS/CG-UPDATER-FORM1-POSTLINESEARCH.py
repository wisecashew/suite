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


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument ("--model", dest='m', type=int, help="Select model directory.")
parser.add_argument ("-T", dest='T', type=float, help="Select temperature.")
args = parser.parse_args()

if __name__=="__main__":

	temps = [args.T]
	for T in temps:
		print (f"In T = {T}", flush=True)
		# get E_mm and delta 
		Emm_old, Ems_old     = aux.get_energy_form_1 (str(T)+"/FORM1/MODEL"+str(args.m)+"/geom_and_esurf.txt")
		delta_cur, scale_old = delta_scale_finder    (str(T)+"/FORM1/MODEL"+str(args.m)+"/delta.mc")
		print (f"delta_cur = {delta_cur}")
		try:
			delta_old, scale_old = delta_scale_finder (str(T)+"/FORM1/MODEL"+str(args.m-1)+"/delta.mc")
			print (f"delta_cur = {delta_cur}, delta_old = {delta_old}")
			if np.sign (delta_cur) * np.sign (delta_old) == -1:
				scale = scale_old*0.5
			else:
				scale = scale_old
		except FileNotFoundError:
			scale = 0.5
			pass
		print ("scale = {:2.10f}".format(scale) )
		Emm_new  = Emm_old - scale * np.sign (delta_cur)
		file_new = open (str(T)+"/FORM1/MODEL"+str(args.m+1)+"/geom_and_esurf.txt", 'w')
		file_new.write ("x = 34\n")
		file_new.write ("y = 34\n")
		file_new.write ("z = 34\n")
		file_new.write (f"kT = {T}\n")
		file_new.write ("Emm = {:2.10f}\n".format (Emm_new) )
		file_new.write ("Ems = 0\n")
		file_new.write ("END OF FILE")
		file_new.close ()
