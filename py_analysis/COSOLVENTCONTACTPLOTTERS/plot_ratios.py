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
parser.add_argument('--contacts'  , dest='contacts',  action='store', type=str, help='enter the address of the real contacts.')
parser.add_argument('--solvation' , dest='solvation', action='store', type=str, help='enter the address of solvation contacts.')
parser.add_argument('--forcefield'  , dest='ff', action='store', type=str, nargs='+', help="enter the forcefield name.")
parser.add_argument('--suffix', dest='s', action='store', type=str, help='enter suffix to images.')
parser.add_argument('--color', dest='color', action='store', type=str, nargs='+', help='enter color of plot.')
parser.add_argument('--show-ylabels', dest='show_ylabels', action='store_true', default=False, help='enter suffix to images.')
parser.add_argument('--scale', dest='scale', action='store_true', default=False, help='enter if you want scaling')
args = parser.parse_args() 

grey_norm  = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.1, vmax=0 ) 
green_norm = matplotlib.colors.Normalize  ( vmin=-3.0, vmax=-0.1 ) 

figsize     = (1.5, 1.5)
marker_size = 8/1.3

def color_finder (H_mix):
	if H_mix <= -0.1:
		return cm.Greens_r (green_norm (H_mix) )
	elif H_mix >-0.1:
		return cm.Greys (grey_norm (H_mix) )

######################################################

def intersection(lst1, lst2):
	return list(set(lst1) & set(lst2))

def plot_contents_of_ss(fidx, df_contacts, df_counts, x_real, colors):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	if not args.show_ylabels:
		ax.set_ylim(0, 0.5)
		ax.set_yticks(np.arange(0, 0.51, 0.1))
		ax.set_yticklabels([])
	frac_s = df_counts["total_s"].values
	frac_c = df_counts["total_c"].values
	ax.plot(x_real, frac_s/(z*M),        marker='o', c=colors[0], linewidth=1,     markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.plot(x_real, frac_c/(z*M),        marker='o', c=colors[1], linewidth=1,     markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.plot(x_real, (frac_s+frac_c)/(z*M), marker='o', c=colors[2], linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.minorticks_on()
	fig.savefig ("solvation-composition-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_ms_a_over_s(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ratio = df_contacts["M1-S1-A"].values/df_counts["total_s"].values
	ratio[-1] = 0
	if not args.show_ylabels:
		ax.set_ylim(0, 3)
		ax.set_yticks(np.arange(0, 3.1, 1))
		ax.set_yticklabels([])
	ax.minorticks_on()
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("aligned-ms-over-total-s-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_ms_a_over_total(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ratio = df_contacts["M1-S1-A"].values/df_counts["total"].values
	if not args.show_ylabels:
		ax.set_ylim(0, 3)
		ax.set_yticks(np.arange(0, 3.1, 1))
		ax.set_yticklabels([])
	ax.minorticks_on()
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("aligned-ms-over-total-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_contacts_over_particles(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ratio = (df_contacts["M1-S1"].values + df_contacts["M1-S2"].values)/(df_counts["total_s"].values+df_counts["total_c"].values)
	ax.minorticks_on()
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("total-contacts-over-total-particles-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_aligned_ms_over_total_ms(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ratio = (df_contacts["M1-S1-A"].values)/(df_contacts["M1-S1"].values)
	ratio[-1] = 0
	ax.minorticks_on()
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("ms1-a-over-total-ms1-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return


def plot_sc_aligned_contacts_over_particles(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ratio = (df_counts["aligned_sc"])/((df_counts["total_s"].values+df_counts["total_c"].values))
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	if not args.show_ylabels:
		ax.set_ylim(0, 10)
		ax.set_yticks(np.arange(0, 10.1, 2))
		ax.set_yticklabels([])
	ax.minorticks_on()
	fig.savefig ("total-sc-aligned-contacts-over-total-particles-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_sc_aligned_contacts_over_total_sc(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	if not args.show_ylabels:
		ax.set_ylim(0, 10)
		ax.set_yticks(np.arange(0, 10.1, 2))
		ax.set_yticklabels([])
	ratio = (df_counts["aligned_sc"].values)/((df_counts["aligned_sc"].values+df_counts["misaligned_sc"].values))
	ratio[0] = 0; ratio[-1] = 0;
	ax.minorticks_on()
	ax.plot (x_real, ratio, marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("total-sc-aligned-contacts-over-total-sc-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_alignment_stats_solvation(fidx, df_contacts, df_counts, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=figsize)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	if not args.show_ylabels:
		ax.set_ylim(0, 0.8)
		ax.set_yticks(np.arange(0, 0.9, 0.2))
		ax.set_yticklabels([])
	ratio_sc = (df_counts["aligned_sc"])/(z*(df_counts["total_s"].values+df_counts["total_c"].values))
	ratio_ms = df_contacts["M1-S1-A"].values/(z*(M)) # df_counts["total_s"].values
	ratio_ms[-1] = 0
	summation = (df_counts["aligned_sc"].values+df_contacts["M1-S1-A"].values)/(z*(df_counts["total_s"].values+df_counts["total_c"].values+M))
	ax.minorticks_on()
	ax.plot (x_real, ratio_sc,   marker='o', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.plot (x_real, ratio_ms,   marker='^', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.plot (x_real, summation,  marker='s', c=color, linewidth=1, markersize=marker_size, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	ax.axvline(x=0.4, color='darkred', linewidth=0.5, ls='--')
	fig.savefig ("total-sc-aligned-contacts-over-total-particles-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

if __name__=="__main__":

	z       = 26
	M       = 32
	ff      = args.ff
	cols    = args.color

	# set up the contacts dataframe
	df_contacts = pd.read_csv(args.contacts, sep='\|', names=["U", "x", "M1-M1", "M1-M1-A", "M1-M1-N", "M1-S", "M1-S1", "M1-S1-A", "M1-S1-N", "M1-S2", "M1-S2-A", "M1-S2-N", "S1-S2", "S1-S2-A", "S1-S2-N"], engine='python', skiprows=1)
	df_contacts = df_contacts[df_contacts["U"] == ff[0]]
	x_real  = df_contacts["x"].values

	# set up the count dataframe
	keys    = ["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "total_sc"]
	df_counts = pd.read_csv(args.solvation, sep='\|', names=["U", "x", "total", "total_s", "total_c", "aligned_sc", "misaligned_sc"], engine='python', skiprows=1)
	df_counts = df_counts[df_counts["U"] == ff[0]]

	# fig, ax = plt.subplots(1, 1, num=1, figsize=(1.5,1.5), squeeze=False)
	plot_contents_of_ss(0, df_contacts, df_counts, x_real, cols)
	plot_ms_a_over_s(0, df_contacts, df_counts, x_real, cols[-1])
	plot_ms_a_over_total(1, df_contacts, df_counts, x_real, cols[-1])
	# plot_contacts_over_particles(1, df_contacts, df_counts, x_real, cols[0])
	plot_sc_aligned_contacts_over_particles(2, df_contacts, df_counts, x_real, cols[-1])
	plot_sc_aligned_contacts_over_total_sc(3, df_contacts, df_counts, x_real, cols[-1])
	plot_aligned_ms_over_total_ms(4, df_contacts, df_counts, x_real, cols[-1])
	plot_alignment_stats_solvation(5, df_contacts, df_counts, x_real, cols[-1])



