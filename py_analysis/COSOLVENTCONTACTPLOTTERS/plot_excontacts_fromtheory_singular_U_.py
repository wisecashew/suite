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
parser.add_argument('--show-ylabels', dest='show_ylabels', action='store_true', default=False, help='enter suffix to images.')
parser.add_argument('--scale', dest='scale', action='store_true', default=False, help='enter if you want scaling')
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

def plot_mm(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ax.set_ylim(-1, 3)
	ax.set_yticks(np.arange(-1, 3.1, 1))
	ax.set_yticklabels([])
	ax.minorticks_on()
	shift     = 63.05 # calculated from CONTACTS-ATHERMAL FOR N=32!!
	df_subset = df_real["M1-M1"].values

	# ax.set_xlabel("$x_c$")
	# ax.set_ylabel("$N_{MM} ^{ex}$")
	ax.plot (x_real, (df_subset-shift)/(M), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-M1-M1-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_ms(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_xlim(0.0, 1.0)
	ax.set_xticks([0, 0.5, 1])
	ax.set_xticklabels([])
	ax.set_ylim(-0.15, 0.1)
	ax.set_yticks(np.arange(-0.15, 0.11, 0.05))
	# ax.set_yticklabels([-0.15, -0.1, -0.05, 0, 0.05, 0.1])
	ax.set_yticklabels([])
	ax.minorticks_on()
	df_subset = df_real ["M1-S1"].values + df_real ["M1-S2"].values
	shift     = (832 - 2 * 63.05)
	df_subset = np.array([val-shift for val in df_subset])
	# ax.set_xlabel("$x_c$")
	# ax.set_ylabel("$N_{MS} ^{ex} + N_{MC} ^{ex}$")
	ax.plot(x_real, (df_subset)/(M*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig("contact-plots-class-M1-S-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_s1s2(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	df_subset = df_real ["S1-S2"].values
	
	# ax.set_xlabel("$x_c$")
	# ax.set_ylabel("$N_{SC} ^{ex}$")
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(0, 0.4)
		ax.set_yticks(np.arange(-0.1, 0.41, 0.1))
		ax.set_yticklabels([])
		ax.minorticks_on()
		contacts = []
		for i in range(len(df_subset)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			contacts.append (df_subset[i] - Ns*z*x_real[i])
		contacts = np.array(contacts)
		
		for i in range(len(contacts)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			if contacts[i] == 0 or x_real[i]==0 or Ns == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(Ns*z*x_real[i])
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_real["S1-S2"].values, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')

	fig.savefig ("contact-plots-class-S1-S2-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_m1s1(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["M1-S1"].values
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(-1, 3)
		ax.set_yticks(np.arange(-1, 3.1, 1))
		ax.set_yticklabels([])
		contacts = []
		for i in range(len(df_subset)):
			contacts.append (df_subset[i] - 705.9 * (1-x_real[i]))
		contacts = np.array(contacts)
		for i in range(len(contacts)):
			if contacts[i] == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(705.9*(1-x_real[i]))
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(0, 1)
		ax.set_yticks(np.arange(0, 1.1, 0.25))
		ax.set_yticklabels([])
		contacts = []
		ax.plot (x_real, df_subset/(M*z), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')

	fig.savefig ("contact-plots-class-M1-S1-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_m1s2(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["M1-S2"].values
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(-1, 3)
		ax.set_yticks(np.arange(-1, 3.1, 1))
		# ax.set_yticklabels([-1, 0, 1, 2, 3])
		ax.set_yticklabels([])
		contacts = []
		for i in range(len(df_subset)):
			contacts.append (df_subset[i] - 705.9 * (x_real[i]))
		contacts = np.array(contacts)
		for i in range(len(contacts)):
			if contacts[i] == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(705.9 * (x_real[i]))
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-M1-S2-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_m1s1_a(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["M1-S1-A"].values
	if args.scale:
		ax.set_xlim(0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(-1, 3)
		ax.set_yticks(np.arange(-1, 3.1, 1))
		ax.set_yticklabels([])
		contacts  = []
		for i in range(len(df_subset)):
			contacts.append (df_subset[i] - 705.9 * (1-x_real[i]))
		contacts = np.array(contacts)
		for i in range(len(contacts)):
			if contacts[i] == 0 or x_real[i] == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(705.9*(1-x_real[i]))
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.set_xlim(0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.plot (x_real, df_subset/(z*M), marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-M1-S1-A-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_m1s1_n(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["M1-S1-N"].values
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(-1, 3)
		ax.set_yticks(np.arange(-1, 3.1, 1))
		ax.set_yticklabels([])
		contacts  = []
		for i in range(len(df_subset)):
			contacts.append (df_subset[i] - 705.9 * (1-x_real[i]))
		contacts = np.array(contacts)
		for i in range(len(contacts)):
			if contacts[i] == 0 or x_real[i] == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(705.9*(1-x_real[i]))
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-M1-S1-N-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_s1s2_a(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["S1-S2-A"].values
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(0, 0.4)
		ax.set_yticks(np.arange(-0.1, 0.41, 0.1))
		ax.set_yticklabels([-0.1, 0, 0.1, 0.2, 0.3, 0.4])
		contacts = []
		for i in range(len(df_subset)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			contacts.append (df_subset[i] - Ns*z*x_real[i])
		contacts = np.array(contacts)
		
		for i in range(len(contacts)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			if contacts[i] == 0 or x_real[i]==0 or Ns == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(Ns*z*x_real[i])
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-S1-S2-A-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return

def plot_s1s2_n(fidx, df_real, x_real, color):
	fig, ax = plt.subplots(1, 1, num=fidx, figsize=(1.5,1.5))
	ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.minorticks_on()
	df_subset = df_real ["S1-S2-N"].values
	if args.scale:
		ax.set_xlim (0.0, 1.0)
		ax.set_xticks([0, 0.5, 1])
		ax.set_xticklabels([])
		ax.set_ylim(0, 0.4)
		ax.set_yticks(np.arange(-0.1, 0.41, 0.1))
		ax.set_yticklabels([-0.1, 0, 0.1, 0.2, 0.3, 0.4])
		contacts = []
		for i in range(len(df_subset)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			contacts.append (df_subset[i] - Ns*z*x_real[i])
		contacts = np.array(contacts)
		
		for i in range(len(contacts)):
			Ns = int(((34**3) - 32)*(1-x_real[i]))
			if contacts[i] == 0 or x_real[i]==0 or Ns == 0:
				contacts[i] = 0
			else:
				contacts[i] = contacts[i]/(Ns*z*x_real[i])
		ax.plot (x_real, contacts, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	else:
		ax.plot (x_real, df_subset, marker='o', c=color, linewidth=1, markersize=8/1.3, markeredgecolor='k', clip_on=False, zorder=10, ls='--')
	fig.savefig ("contact-plots-class-S1-S2-N-"+args.s, bbox_inches='tight', dpi=1200)
	ax.cla()
	return


if __name__=="__main__":

	z       = 26
	M       = 32
	ff    = args.ff
	cols    = args.color
	keys    = ["M1-M1", "M1-S", "S1-S2", "M1-S1", "M1-S2", "M1-S1-A", "M1-S1-N", "S1-S2-A", "S1-S2-N"]

	# df_ath  = pd.read_csv("CONTACTS-ATHERMAL.csv", sep='\|', names=["H", "x", "M1-M1", "M1-S1", "M1-S2", "S1-S1", "S1-S2", "S2-S2"], engine='python', skiprows=1)
	df_real = pd.read_csv(args.real, sep='\|', names=["H", "x", "M1-M1", "M1-M1-A", "M1-M1-N", "M1-S", "M1-S1", "M1-S1-A", "M1-S1-N", "M1-S2", "M1-S2-A", "M1-S2-N", "S1-S2", "S1-S2-A", "S1-S2-N"], engine='python', skiprows=1)
	df_real = df_real[df_real["H"] == ff[0]]
	x_real  = df_real["x"].values

	# fig, ax = plt.subplots(1, 1, num=1, figsize=(1.5,1.5), squeeze=False)
	for idx,k in enumerate(keys):
		if k == "M1-M1":
			plot_mm(idx, df_real, x_real, color=args.color[0])
		elif k == "M1-S":
			plot_ms(idx, df_real, x_real, color=args.color[0])
		elif k == "S1-S2":
			plot_s1s2(idx, df_real, x_real, color=args.color[0])
		elif k == "M1-S1":
			plot_m1s1(idx, df_real, x_real, color=args.color[0])
		elif k == "M1-S2":
			plot_m1s2(idx, df_real, x_real, color=args.color[0])
		elif k == "M1-S1-A":
			plot_m1s1_a(idx, df_real, x_real, color=args.color[0])
		elif k == "M1-S1-N":
			plot_m1s1_n(idx, df_real, x_real, color=args.color[0])
		elif k == "S1-S2-A":
			plot_s1s2_a(idx, df_real, x_real, color=args.color[0])
		elif k == "S1-S2-N":
			plot_s1s2_n(idx, df_real, x_real, color=args.color[0])
		else:
			print(f"Invalid key: {k}! Exiting...")
			exit()



