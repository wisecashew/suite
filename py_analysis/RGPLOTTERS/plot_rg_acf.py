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
parser.add_argument('--nproc',        metavar='N',    type=int,              dest='nproc',   action='store', help='Request these many proccesses.')
parser.add_argument('--frac',         metavar='frac', dest='frac',           action='store', nargs='+',      type=float, help='Enter fractions you want probed.')
parser.add_argument('--U',            metavar='U',    dest='U',              action='store', nargs='+',      type=str,   help='Enter U you want probed.')
parser.add_argument('--energy',       metavar='dump', type=str,              dest='dump',    action='store', help='Energy file to look for starting index.')
parser.add_argument('--coords',       dest='c',       metavar='coords',      action='store', type=str,       help='Name of energy dump file to parse information.', default='coords.txt')
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
	print(cols, flush=True)
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

	nproc = args.nproc
	pool = multiprocessing.Pool ( processes=nproc )# len(num_list)) 

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
			num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(f) ), args.dump ) ) )
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
		idx_range = len (master_num_list)//nproc + 1
		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				results = pool.starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_frac_list[u_idx*nproc:], master_num_list [u_idx*nproc:], itertools.repeat (dop), itertools.repeat (coords_files), master_index_list[u_idx*nproc:] ) )
			else:
				results = pool.starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_frac_list[(u_idx)*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(coords_files), master_index_list[u_idx*nproc:(u_idx+1)*nproc] ) )
			print ("\tPool has been closed. This pool had {} threads.".format (len(results) ), flush=True )
			for k in range( len( master_frac_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				rg_dict[master_frac_list[u_idx*nproc + k]].append( results[k] )

		# collect the rgs
		acf_list = []
		for idx,f in enumerate(frac_list):
			print(f"@ fraction = {f}", flush=True)
			print(f"This should be empty: {acf_list}", flush=True)
			print(f"Length of rg_dict[{f}] = {len(rg_dict[f])}", flush=True)
			for rg_list in rg_dict[f]:
				autocorr = acf(rg_list, nlags=100, fft=True)
				acf_list.append(autocorr)
			acfs = np.array(acf_list)
			acf_mean = np.mean(acfs, axis=0)
			ax.plot(range(101), acf_mean, marker='o', mec='k', ls='--', c=cols[idx], linewidth=1, ms=8/1.3, clip_on=False, label="$x_{\\mathrm{c}}$ = " + str(f))
			acf_list.clear()

	# close the pool
	pool.close()
	pool.join()

	# set up the plots
	ax.set_ylim(-0.1,1)
	ax.set_xlim(0, 100)
	ax.legend(loc='upper right', fontsize=8, frameon=False)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	print(f"Made autocorrelation.", flush=True)
	fig.savefig(args.pn+".png", dpi=1200, bbox_inches="tight")

	##################################
	stop = time.time() 
	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)
