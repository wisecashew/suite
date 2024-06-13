#!/home/satyend/.conda/envs/phase/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis")
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
parser.add_argument('--real'  , dest='real', action='store', type=str, help='enter the address of the real contacts.')
parser.add_argument('--colors', dest='cols', action='store', nargs='+', type=str, help='Enter the colors you want.')
parser.add_argument('--U'  , dest='U', action='store', type=str, nargs='+', help="enter the enthalpic conditions.")
parser.add_argument('--suffix', dest='s', action='store', type=str, help='enter suffix to images.')
parser.add_argument('--show-ylabels', dest='show_ylabels', action='store_true', default=False, help='enter suffix to images.')
args = parser.parse_args() 

grey_norm  = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.1, vmax=0 ) 
green_norm = matplotlib.colors.Normalize  ( vmin=-3.0, vmax=-0.1 ) 

######################################################

def intersection(lst1, lst2):
	return list(set(lst1) & set(lst2))

if __name__=="__main__":

	# cols = ["rosybrown", "lightcoral", "indianred", "brown"]
	z       = 26
	M       = 32
	U_list  = args.U
	cols    = args.cols
	fig, ax = plt.subplots (1, 1, num=1, squeeze=False, figsize=(2.5,2.5))

	keys    = ["SOLVATION"]

	df_real = pd.read_csv (args.real, sep='|', names=["U", "x", "SOLVATION", "T_ms", "N_ms"], engine='python', skiprows=1)
	x_real  = df_real ["x"]

	ax[0][0].tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax[0][0].tick_params(axis='x', labelsize=8)
	ax[0][0].tick_params(axis='y', labelsize=8)
	if args.show_ylabels:
		pass
	else:
		ax[0][0].set_yticklabels([])
	ax[0][0].set_xlim (0.0, 1.0)
	ax[0][0].set_xticks (np.arange(0,1.2,0.2))
	ax[0][0].set_xticklabels ([]) # ax[0][0].get_xticks(), weight='bold')
	ax[0][0].minorticks_on()

	for idx, hmix in enumerate(U_list):
		df_subset = df_real [df_real ["U"] == hmix]
		idx       = idx % len(cols)
		ax[0][0].plot (df_subset["x"].values, df_subset["SOLVATION"].values, marker='o', c=cols[idx], linewidth=1, markersize=8/1.3, markeredgecolor='k', label=f"{hmix}", clip_on=False, zorder=10)

	plt.savefig ("solvation"+"-"+args.s, bbox_inches='tight', dpi=1200)
