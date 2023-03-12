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

	fig = plt.figure( figsize=(4/2,3/2), constrained_layout=True )
	ax  = plt.axes() 
	plt.rcParams["axes.labelweight"] = "bold"
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	
	temps = aux.dir2float (os.listdir ("."))
	temps.remove (0.09)
	temps.remove (0.3)
	chi_list = np.zeros (len(temps))
	count = 0
	for T in temps:
		print (f"In T = {T}", flush=True)
		# get E_mm and delta 
		Emm, Ems = aux.get_energy_cg (str(T)+"/MODEL"+str(args.m)+"/geom_and_esurf.txt")
		chi_list[count] = (Ems - 0.5*Emm)/T
		count += 1
	ax.set_xscale ("log")
	ax.minorticks_on()
	ax.set_xticks ([0.01, 0.1, 1, 10, 100])
	ax.set_yticklabels (ax.get_yticks(), weight='bold')
	ax.set_xticklabels ([0.01, 0.1, 1, 10, 100], weight='bold')
	plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
	plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
	plt.plot (temps, chi_list, marker='o', markeredgecolor='k')
	plt.savefig ("chi_over_T", dpi=1200)
