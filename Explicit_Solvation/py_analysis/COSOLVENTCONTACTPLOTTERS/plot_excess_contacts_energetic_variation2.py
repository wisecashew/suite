#!/usr/licensed/anaconda3/2020.7/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import sys
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis')
import aux 
import pandas as pd
import argparse
import itertools
import multiprocessing
import os
import matplotlib.cm as cm

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

parser = argparse.ArgumentParser(description="Go into lattice dump and get contacts.")
parser.add_argument('--ideal', dest='id', type=str, action='store', help='enter the address of the ideal contacts.')
parser.add_argument('--real', dest='real', action='store', type=str, help='enter the address of the real contacts.')
args = parser.parse_args() 

divnorm = matplotlib.colors.Normalize ( vmin=-3, vmax=0 ) 

######################################################

def intersection(lst1, lst2):
	return list(set(lst1) & set(lst2))

if __name__=="__main__":
	df_id   = pd.read_csv (args.id, sep='|')
	x_id    = df_id ["x"]
	plt.rcParams["axes.labelweight"] = "bold"
	keys    = ["M1-M1", "S1-S2", "M1-S1", "M1-S2", "S1-S1", "S2-S2"]
	titles  = ["Monomer-monomer contacts", "Solvent-cosolvent contacts", "Monomer-solvent contacts", "Monomer-cosolvent contacts", "Solvent-solvent contacts", "cosolvent-cosolvent contacts"]
	ylims   = [(-0.55,2.1),(-0.001,0.35), (-1.125, 0.375), (-1.125, 0.375), (-1.125, 0.375), (-1.125, 0.375)]
	ypads   = [3, 1, 1, 1, 1, 3]
	yticks  = [np.arange(-0.5, 2.5, 0.5), np.arange(0, 0.4, 0.05), np.arange(-1.0, 0.375, 0.25), np.arange(-1.0, 0.375, 0.25), np.arange(-1.0, 0.5, 0.25), np.arange (-1, 0.375, 0.25)]
	fig, ax = plt.subplots (1, 1, num=1, squeeze=False, figsize=(4/2,3/2))
	
	df_real = pd.read_csv (args.real, sep='|')
	x_real  = df_real ["x"]
	frac = intersection (x_id, x_real)
	df_id = df_id [df_id["x"].isin (frac)]
	enthalpies = np.unique(df_real["H"])
	print (df_id)
	colormap = cm.get_cmap("coolwarm")
	for k in range( 4 ): # len(keys)):
		for i in range(len(enthalpies)):
			df_real_subset = df_real [df_real ["H"] == enthalpies[i]]
			df_real_subset = df_real_subset[df_real_subset["x"].isin (frac)]
			print (df_real_subset)
			ax[0][0].tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
			ax[0][0].tick_params(axis='x', labelsize=8) # , pad=3)
			ax[0][0].tick_params(axis='y', labelsize=8) # , pad=ypads[k])
			# ax[0][0].set_ylim (ylims[k])
			# ax[0][0].set_yticks (yticks[k])
			ax[0][0].set_yticklabels (ax[0][0].get_yticks(), weight='bold')
			ax[0][0].set_xlim (-0.05, 1.05)
			ax[0][0].set_xticks (np.arange(0,1.2,0.2))
			ax[0][0].set_xticklabels (ax[0][0].get_xticks(), weight='bold')
			# ax[0][0].set_xlabel ("$x_c$", fontsize=8, labelpad=2)
			ax[0][0].minorticks_on()
			if k == 0 or k > 1:
				ax[0][0].yaxis.set_major_formatter (matplotlib.ticker.StrMethodFormatter ('{x:1.1f}'))
			elif k == 1:
				ax[0][0].yaxis.set_major_formatter (matplotlib.ticker.StrMethodFormatter ('{x:1.2f}'))
			ax[0][0].xaxis.set_major_formatter (matplotlib.ticker.StrMethodFormatter ('{x:1.1f}'))
			p_diff = (df_real_subset[keys[k]].values - df_id[keys[k]].values)/(df_real_subset[keys[k]].values)
			ax[0][0].plot (df_id["x"].values, p_diff, marker='o', linewidth=3/2, markersize=8/2, color=cm.summer(divnorm(enthalpies[i])), markeredgecolor='k')
	
		plt.savefig ("excess-plots-class-"+keys[k], bbox_inches='tight', dpi=1200)
		ax[0][0].cla()
