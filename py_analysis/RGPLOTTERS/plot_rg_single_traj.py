#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import os
import sys
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux
import time
import multiprocessing
from statsmodels.tsa.stattools import acf
import itertools

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

'''
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

import argparse
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of fraction and potential energy surfaces.")
parser.add_argument('-dop',           metavar='DOP',   dest='dop',            type=int,      action='store', help='enter a degree of polymerization.')
parser.add_argument('-s',             metavar='S',     type=int,              dest='s',      action='store', help='start parsing after this move number (not index or line number in file).', default=100)
parser.add_argument('--frac',         metavar='frac', dest='frac',           action='store', nargs='+',      type=float, help='Enter fractions you want probed.')
parser.add_argument('--U',            metavar='U',    dest='U',              action='store', nargs='+',      type=str,   help='Enter U you want probed.')
parser.add_argument('--energy',       metavar='dump', type=str,              dest='dump',    action='store', help='Energy file to look for starting index.')
parser.add_argument('--coords',       dest='c',       metavar='coords',      action='store', type=str,       help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--cnum',         dest='cnum',    metavar='coords',      action='store', type=int,       help='Number of traj to probe.', default='coords.txt')
parser.add_argument('--mean',         metavar='mean', dest='mean',           action='store', type=float, help='Enter mean.')
parser.add_argument('--figsize',      dest='fs',      metavar='figsize',   action='store', type=float,  nargs='+', help='Size of figure.', default=[2.5, 2.5])
parser.add_argument('--colors',  dest='colors', action='store', nargs='+', type=str,  help='Color of plot.', default=False)
parser.add_argument('--png-name',     dest='pn',      metavar='imagename',   action='store', type=str,       help='Name of image file', default='rg_plot')
args = parser.parse_args() 

def get_starting_ind ( U, frac, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + str(frac) + "/" + dumpfile + "_" + str(num) + ".mc"
	df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	L = len(df["energy"])
	return int(df["time_step"].values[L-args.s])


def infiltrate_coords_get_rg ( U, frac, num, dop, coords_files, starting_index ):

	filename = U + "/DOP_" + str(dop) + "/" + str(frac) + "/"+ coords_files + "_" + str(num)+".mc" 
	edge = aux.edge_length (dop)
	master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
	rg = aux.get_Rg_list(master_dict, edge, edge, edge) 
	# print(f"Rg mean = {np.mean(rg)}", flush=True)
	return rg 


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

def mysorter_f (x):
	new_str = ""
	for m in x:
		if m.isdigit():
			new_str = new_str+m
	return float(new_str)

if __name__ == "__main__":


	start = time.time()
	U_list = args.U
	cols    = args.colors # color_map("coral", "darkred", len(U_list))
	print (U_list, flush=True)
	PLOT_DICT = {} 
	dop            = args.dop
	coords_files   = args.c
	starting_index = args.s

	fig = plt.figure( figsize=(args.fs[0],args.fs[1]))
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	i = 0 

	for U in U_list:
		print (f"Diving into U = {U}...", flush=True)
		frac_list = args.frac
		rg_mean = [] 
		rg_std  = [] 
		
		# get num_list for each temperature 
		master_frac_list     = []
		master_num_list      = []
		master_index_list    = []
		rg_dict    = {}
		ntraj_dict = {}
		for f in frac_list:
			num_list = [args.cnum]
			master_num_list.extend ( num_list )
			for num in num_list:
				master_index_list.append (get_starting_ind (U, f, num, dop, args.dump) )

			master_frac_list.extend ( [f]*len( num_list ) )
			ntraj_dict[f] = len ( num_list )
			rg_dict[f]    = []

		if len(num_list) == 0:
			print(f"There is an issue with dir2nsim. Check it out.", flush=True)
			exit()

		# start multiprocessing... keeping in mind that each node only has 96 cores 
		# start splitting up master_num_list and master_temp_list 
		idx_range = len (master_num_list)//1 + 1
		my_rg_list = infiltrate_coords_get_rg(U, master_frac_list[0], master_num_list[0], dop, coords_files, master_index_list[0])
		my_rg_list = np.array(my_rg_list)/args.mean
	for idx,U in enumerate(U_list):
		ax.plot(range(len(my_rg_list)), my_rg_list, marker='o', markeredgecolor='k', linestyle='--', \
			c=cols[0], linewidth=1, markersize=8/1.3, label=f'{U}', clip_on=False, zorder=10)
		i += 1
		ax.minorticks_on()
		ax.set_ylim(0, 3)
		ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
		ax.set_xlim(0, len(my_rg_list)-1)
		ax.set_yticklabels([])
		ax.set_xticklabels([])
		print("Made Rg timeseries.", flush=True)
		plt.savefig(args.pn+"-timeseries.png", bbox_inches='tight', dpi=1200)

	n_blocks   = 8
	block_size = len(my_rg_list) // n_blocks
	block_means = np.array([np.mean(my_rg_list[i*block_size:(i+1)*block_size]) for i in range(n_blocks)])
	block_std   = np.array([np.std (my_rg_list[i*block_size:(i+1)*block_size]) for i in range(n_blocks)])
	fig = plt.figure(figsize=(args.fs[0], args.fs[1]))
	ax  = plt.axes()
	ax.tick_params(direction="in", top=True, right=True, bottom=True, left=True)
	ax.axhline(y=np.mean(block_means), c="lavender", ls='--')
	ax.errorbar(range(n_blocks), block_means, yerr=block_std, marker='^', ecolor='k', markeredgecolor='k', linestyle='-', c=cols[1], linewidth=1, ms=8/1.3, clip_on=False)
	print(f"Obtained block averaged.", flush=True)
	ax.set_ylim(0, 2)
	ax.set_yticks([0, 0.5, 1, 1.5, 2])
	ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7])
	ax.set_xlim(0, n_blocks-1)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	fig.savefig(args.pn+"-blockave.png", dpi=1200, bbox_inches="tight")
	fig = plt.figure(figsize=(args.fs[0], args.fs[1]))
	ax  = plt.axes()
	ax.tick_params(direction="in", top=True, right=True, bottom=True, left=True)
	autocorr = acf(my_rg_list, nlags=100, fft=True)
	ax.plot(range(101), autocorr, marker='o', mec='k', ls='--', c=cols[3], linewidth=1, ms=8/1.3, clip_on=False)
	ax.set_ylim(-0.1,1)
	ax.set_xlim(0, 100)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	print(f"Made autocorrelation.", flush=True)
	fig.savefig(args.pn+"-acf.png", dpi=1200, bbox_inches="tight")

	##################################
	stop = time.time() 
	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)
