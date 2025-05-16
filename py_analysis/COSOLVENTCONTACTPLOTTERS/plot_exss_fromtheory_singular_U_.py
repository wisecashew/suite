#!/home/satyend/.conda/envs/phase/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use("agg")
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
parser.add_argument('--forcefield'  , dest='ff', action='store', type=str, nargs='+', help="enter the forcefield name.")
parser.add_argument('--suffix', dest='s', action='store', type=str, help='enter suffix to images.')
parser.add_argument('--color', dest='color', action='store', type=str, nargs='+', help='enter color of plot.')
parser.add_argument('--scale', dest='scale', action='store_true', default=False, help='enter if you want scaling')
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

def plot_total(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	# ax.set_xticklabels([])
	# ax.set_ylim(-1, 3)
	# ax.set_yticks(np.arange(-1, 3.1, 1))
	# ax.set_yticklabels([])
	ax.minorticks_on()
	df_subset = df_real["total"].values
	if args.scale:
		ax.plot (x_real, df_subset/(M*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("number-plots-total-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_total_s(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	# ax.set_ylim(-0.15, 0.1)
	# ax.set_yticks(np.arange(-0.15, 0.11, 0.05))
	# ax.set_yticklabels([-0.15, -0.1, -0.05, 0, 0.05, 0.1])
	# ax.set_yticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["total_s"].values
	if args.scale:
		to_plot = [val/(M*z*(1-x)) if x != 1 else 0 for val,x in zip(df_subset, x_real)]
		ax.plot (x_real, to_plot, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.set_ylim(0, 0.5)
		ax.set_yticks(np.arange(0.0, 0.51, 0.1))
		ax.set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5])
		ax.set_yticklabels([])
		ax.plot (x_real, df_subset/(M*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig("number-plots-total-s-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_total_c(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim (0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	# ax.set_ylim(0, 0.4)
	# ax.set_yticks(np.arange(-0.1, 0.41, 0.1))
	# ax.set_yticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["total_c"].values
	if args.scale:
		to_plot = [val/(M*z*x) if x != 0 else 0 for val,x in zip(df_subset, x_real)]
		ax.plot (x_real, to_plot, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.set_ylim(0, 0.5)
		ax.set_yticks(np.arange(0.0, 0.51, 0.1))
		ax.set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5])
		ax.set_yticklabels([])
		ax.plot (x_real, df_subset/(M*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("number-plots-total-c-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_total_sc(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim (0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["aligned_sc"].values +df_real ["misaligned_sc"].values
	if args.scale:
		ax.set_ylim(0, 0.5)
		ax.set_yticks(np.arange(0, 0.55, 0.1))
		ax.set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5])
		ax.set_yticklabels([])
		to_plot = [val/(M*z*x*z*(1-x)) if (x != 0 and x != 1) else 0 for val,x in zip(df_subset, x_real)]
		ax.plot (x_real, to_plot, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("number-plots-total-SC-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_aligned_sc(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim (0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["aligned_sc"].values
	if args.scale:
		ax.set_ylim(0, 0.5)
		ax.set_yticks(np.arange(0, 0.55, 0.1))
		ax.set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5])
		ax.set_yticklabels([])
		to_plot = [val/(M*z*x*z*(1-x)) if (x != 0 and x != 1) else 0 for val,x in zip(df_subset, x_real)]
		ax.plot (x_real, to_plot, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.set_ylim(0, 0.2)
		ax.set_yticks(np.arange(0, 0.21, 0.05))
		ax.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
		ax.set_yticklabels([])
		ax.plot (x_real, df_subset/(M*z*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("number-plots-aligned-SC-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_misaligned_sc(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim (0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	# ax.set_xticklabels([])
	# ax.set_ylim(-1, 3)
	# ax.set_yticks(np.arange(-1, 3.1, 1))
	# ax.set_yticklabels([-1, 0, 1, 2, 3])
	# ax.set_yticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["misaligned_sc"].values
	if args.scale:
		to_plot = [val/(M*z*x*z*(1-x)) if (x != 0 and x != 1) else 0 for val,x in zip(df_subset, x_real)]
		ax.plot (x_real, to_plot, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("number-plots-misaligned-SC-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

if __name__=="__main__":

	z       = 26
	M       = 32
	ff    = args.ff
	cols    = args.color
	keys    = ["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "total_sc"]

	df_real = pd.read_csv(args.real, sep='\|', names=["U", "x", "total", "total_s", "total_c", "aligned_sc", "misaligned_sc"], engine='python', skiprows=1)
	df_real = df_real[df_real["U"] == ff[0]]
	x_real  = df_real["x"].values

	# fig, ax = plt.subplots(1, 1, num=1, figsize=(1.5,1.5), squeeze=False)
	for idx,k in enumerate(keys):
		if k == "total":
			plot_total(idx, df_real, x_real, color=args.color[0])
		elif k == "total_s":
			plot_total_s(idx, df_real, x_real, color=args.color[0])
		elif k == "total_c":
			plot_total_c(idx, df_real, x_real, color=args.color[0])
		elif k == "aligned_sc":
			plot_aligned_sc(idx, df_real, x_real, color=args.color[0])
		elif k == "misaligned_sc":
			plot_misaligned_sc(idx, df_real, x_real, color=args.color[0])
		elif k == "total_sc":
			plot_total_sc(idx, df_real, x_real, color=args.color[0])
		else:
			print(f"Invalid key: {k}! Exiting...")
			exit()

