#!/usr/licensed/anaconda3/2020.7/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import aux 
import pandas as pd
import argparse
import itertools
import multiprocessing
import os
import sys
import matplotlib.cm as cm

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

parser = argparse.ArgumentParser(description="Go into lattice dump and get contacts.")
parser.add_argument('--ideal', dest='id', type=str, action='store', help='enter the address of the ideal contacts.')
parser.add_argument('--real', dest='real', action='store', type=str, help='enter the address of the real contacts.')
parser.add_argument('-T', dest='T', action='store', nargs='+', type=float, help='enter temperatures at which we are going to do comparisons.')
args = parser.parse_args() 

######################################################
if __name__=="__main__":
	df_id = pd.read_csv (args.id, sep='|')
	keys = ["M1-M1", "S1-S2", "M1-S1", "M1-S2", "S1-S1", "S2-S2"]
	titles = ["Monomer-monomer contacts", "Solvent-cosolvent contacts", "Monomer-solvent contacts", "Monomer-cosolvent contacts", "Solvent-solvent contacts", "Solvent-cosolvent contacts"]
	ylims  = [(-1.2,2.2), (-0.1,0.4), (-1.5, 0.7), (-1.5, 0.7), (-1.25, 0.5), (-1.25, 0.5)]
	ypads  = [3, 1, 1, 1, 1, 3]
	yticks = [np.arange(-1, 2.5, 0.5), np.arange(0, 0.4, 0.05), np.arange(-1.25, 0.5+0.25, 0.25), np.arange(-1.25, 0.75, 0.25), np.arange(-1.0, 0.5, 0.25)]
	# yticks = [np.arange(-1.5, 2.5, 0.5), np.arange(-1.5, 0.50, 0.3), (-1.5, 0.5), (-1.5, 0.5), (-0.2, 0.5), (-1.5, 0.5)]
	fig, ax = plt.subplots (len(args.T), len(keys), squeeze=False, figsize=(10,2*len(args.T)))
	
	Tdict = {}
	colormap = cm.get_cmap("jet")
	for k in range(len(keys)):
		Tdict.clear() 
		for T in args.T:
			Tdict[T] = pd.read_csv (args.real+"_T_"+str(T)+".csv", sep='|')
	
		for i in range(len(args.T)):
			if i == 0:
				ax[i][k].set_title (titles[k], fontsize=6)
			if k == 0:
				ax[i][k].set_ylabel ("% difference", fontsize=5, labelpad=2)
			ax[i][k].tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
			ax[i][k].tick_params(axis='x', labelsize=4, pad=3)
			ax[i][k].tick_params(axis='y', labelsize=4, pad=ypads[i])
			ax[i][k].set_ylim (ylims[k])
			ax[i][k].set_yticks (yticks[k])
			ax[i][k].set_xlim (-0.1, 1.1)
			ax[i][k].set_xticks (np.arange(0,1.2,0.2))
			ax[i][k].set_xlabel ("$x_c$", fontsize=5, labelpad=2)
			if np.mean (df_id[keys[k]].values) == 0.0:
				p_diff = np.zeros (len(df_id[keys[k]].values))
			else:
				p_diff = (Tdict[args.T[i]][keys[k]].values - df_id[keys[k]].values)/(df_id[keys[k]].values)
			ax[i][k].plot (df_id["x"].values, p_diff, marker='o', linewidth=2, markersize=6, color=colormap(k/len(keys)))
		
	plt.savefig ("plots.png", bbox_inches='tight', dpi=1200)

