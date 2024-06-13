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

def color_finder (H_mix):
	if H_mix <= -0.1:
		return cm.Greens_r (green_norm (H_mix) )
	elif H_mix >-0.1:
		return cm.Greys (grey_norm (H_mix) )

######################################################

def intersection(lst1, lst2):
	return list(set(lst1) & set(lst2))

import matplotlib.colors as mcolors
def color_map(start_color, end_color, n_steps):
	start_rgb = mcolors.hex2color(mcolors.CSS4_COLORS[start_color])
	end_rgb   = mcolors.hex2color(mcolors.CSS4_COLORS[end_color])

	# linearly interpolate the RGB values
	r = [start_rgb[0] + (end_rgb[0] - start_rgb[0]) * i/n_steps for i in range(n_steps+1)]
	g = [start_rgb[1] + (end_rgb[1] - start_rgb[1]) * i/n_steps for i in range(n_steps+1)]
	b = [start_rgb[2] + (end_rgb[2] - start_rgb[2]) * i/n_steps for i in range(n_steps+1)]
	
	colors = [mcolors.to_hex([r[i], g[i], b[i]]) for i in range(n_steps+1)]

	return colors


if __name__=="__main__":

	# cols = ["rosybrown", "lightcoral", "indianred", "brown"]
	z       = 26
	M       = 32
	Hmix    = args.U
	cols    = args.cols # color_map("coral", "darkred", len(Hmix))
	fig, ax = plt.subplots (1, 1, num=1, squeeze=False, figsize=(2.5,2.5))

	keys    = ["M1-M1", "M1-S", "S1-S2", "M1-S1", "M1-S2", "M1-S1-A", "M1-S1-N"]
	titles  = ["Monomer-monomer contacts", "Solvent-cosolvent contacts", "Monomer-solvent contacts", "Monomer-cosolvent contacts", "Solvent-solvent contacts", "cosolvent-cosolvent contacts"]

	ylims   = [(-0.55,2.1),(-0.001,0.35), (-1.125, 0.375), (-1.125, 0.375), (-1.125, 0.375), (-1.125, 0.375)]
	ypads   = [3, 1, 1, 1, 1, 3]
	yticks  = [np.arange(-0.5, 2.5, 0.5), np.arange(0, 0.4, 0.05), np.arange(-1.0, 0.375, 0.25), np.arange(-1.0, 0.375, 0.25), np.arange(-1.0, 0.5, 0.25), np.arange (-1, 0.375, 0.25)]

	df_real = pd.read_csv (args.real, sep='|', names=["H", "x", "M1-M1", "M1-M1-A", "M1-M1-N", "M1-S", "M1-S1", "M1-S1-A", "M1-S1-N", "M1-S2", "M1-S2-A", "M1-S2-N", "S1-S2", "S1-S2-A", "S1-S2-N"], engine='python', skiprows=1)
	x_real  = df_real ["x"]

	colormap = cm.get_cmap ("coolwarm")
	for k in range( len(keys) ):

		ax[0][0].tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
		ax[0][0].tick_params(axis='x', labelsize=8)
		ax[0][0].tick_params(axis='y', labelsize=8)
		if args.show_ylabels:
			pass # ax[0][0].get_yticks(), weight='bold')
		else:
			ax[0][0].set_yticklabels([])
		ax[0][0].set_xlim (0.0, 1.0)
		ax[0][0].set_xticks (np.arange(0,1.2,0.2))
		# ax[0][0].set_xticklabels ([]) # ax[0][0].get_xticks(), weight='bold')
		ax[0][0].minorticks_on()
		if k == 0 or k > 1:
			pass 
		elif k == 1:
			pass

		for idx, hmix in enumerate(Hmix):
			df_subset = df_real [df_real ["H"] == hmix]
			print (df_subset)
			idx = idx % len(cols)
			if keys == "M1-S1" or keys == "M1-S2":
				p_diff = (df_subset [keys[k]].values - 0) # df_id ["x"].values*z*M)
			else:
				p_diff = (df_subset [keys[k]].values - 0) # df_id [keys[k]].values)
			ax[0][0].plot (df_subset["x"].values, p_diff, marker='o', c=cols[idx], linewidth=1, markersize=8/1.3, markeredgecolor='k', label=f"{hmix}", clip_on=False, zorder=10)

		plt.savefig ("contact-plots-class-"+keys[k]+"-"+args.s, bbox_inches='tight', dpi=1200)
		ax[0][0].cla()
