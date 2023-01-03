#!/usr/licensed/anaconda3/2020.7/bin/python

import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis")
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse
import aux
import os


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('--database', dest='db', action='store', type=str, help='Name of database.')
parser.add_argument('--name'    , dest='n', action='store', type=str, help="name of image.")
parser.add_argument('--excl-vol', dest='ev', action='store_true', default=False)
args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	U_list = aux.dir2U (os.listdir("."))
	i = 0
	chi_list = [0.1, 0.05, 0.01, 0.005, 0.001, 0, -0.001, -0.01, -0.05, -0.1, -0.2]
	rgba_color = cm.PiYG (divnorm(chi_list[i]))
	fig = plt.figure( num=1, figsize=(8,6) )
	ax  = plt.axes()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
	df = pd.read_csv (args.db, sep='|')
	temperatures = [0.01, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0] # df[ df["U"] == U ]["T"]
	df = df[df["T"].isin (temperatures)]
	for U in U_list:
		rgba_color   = cm.PiYG ( divnorm (chi_list[i]) )
		Rg_mean      = df[ df["U"] == U ]["Rg_mean"]
		Rg_err       = df[ df["U"] == U ]["Rg_err" ]
		ax.errorbar ( temperatures, Rg_mean, yerr=Rg_err, linewidth=1, capsize=2, \
		color=rgba_color, ecolor='k', fmt="none", label="_nolegend_")
		ax.plot ( temperatures, Rg_mean, marker='o', markeredgecolor='k', \
		linestyle='-', linewidth=3, color=rgba_color, label="_nolegend_", markersize=10 )
		i += 1
	
	if args.ev:
		temperatures = df[ df["U"] == "Uexcl" ]["T"]
		Rg_mean      = df[ df["U"] == "Uexcl" ]["Rg_mean"]
		Rg_err       = df[ df["U"] == "Uexcl" ]["Rg_err"]
		ax.errorbar (temperatures, Rg_mean, yerr=0, linestyle='-', linewidth=3)
	
	ax.set_xscale ("log")
	ax.set_ylim   (0.2, 1.6)
	ax.tick_params(axis='x', labelsize=16)
	ax.tick_params(axis='y', labelsize=16)
	ax.set_yticks (np.arange(0.3, 1.6, 0.3))
	plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.2f}'))
	ax.minorticks_on ()
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	plt.savefig (args.n, bbox_inches='tight', dpi=1200)


