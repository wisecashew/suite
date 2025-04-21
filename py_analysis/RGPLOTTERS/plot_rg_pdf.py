#!/home/satyend/.conda/envs/phase/bin/python

import itertools
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import time
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tck
import argparse
import sys
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux
import os
import multiprocessing
from pathlib import Path
from scipy.stats import beta
from scipy.stats import chi2
from scipy.stats import poisson
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-dop',           metavar='DOP',  dest='dop',            type=int,       action='store', help='enter a degree of polymerization.')
parser.add_argument('-s',             metavar='S',    type=int,              dest='s',       action='store', help='start parsing after this move number (not index or line number in file).', default=100)
parser.add_argument('--T',            metavar='T',    dest='T',           action='store', nargs='+',      type=float, help='Enter the temperatures you want probed.')
parser.add_argument('--U',            metavar='U',    dest='U',              action='store', nargs='+',      type=str,   help='Enter U you want probed.')
parser.add_argument('-nproc',         metavar='N',    type=int,              dest='nproc',   action='store', help='Request these many proccesses.')
parser.add_argument('--coords',       dest='c',       metavar='coords.txt',  action='store', type=str,       help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--show-legends', dest='sl',      action='store_true',   help='Name of energy dump file to parse information.', default=False)
parser.add_argument('--png-name',     dest='pn',      metavar='imagename',   action='store', type=str,       help='Name of image file', default='rg_plot')
args = parser.parse_args() 

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

def get_starting_ind ( U, T, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/" + dumpfile + "_" + str(num) + ".mc"
	df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	L = len(df["energy"])
	return int(df["time_step"].values[L-args.s])


def infiltrate_coords_get_rg ( U, T, num, dop, coords_files, starting_index ):

	filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num)+".mc" 
	edge = aux.edge_length (dop)
	master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
	rg = aux.get_Rg(master_dict, edge, edge, edge) 
	return rg 



if __name__ == "__main__":

	start = time.time() 
	##################################
	def mysorter_f (x):
		new_str = ""
		for m in x:
			if m.isdigit():
				new_str = new_str+m
		return float(new_str)

	U_list = args.U
	print (U_list, flush=True)
	PLOT_DICT = {}
	dop            = args.dop
	coords_files   = args.c
	starting_index = args.s

	######
	fig = plt.figure( figsize=(2.5,2.5), constrained_layout=True )
	ax  = plt.axes() 
	plt.rcParams["axes.labelweight"] = "bold"
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	temperature_list = args.T
	i = 0

	rg_max = (dop ** (0.588))/2 # (1+np.sqrt(2)+np.sqrt(3))/(3*6**0.5) * (dop**0.57) 
	# instantiating pool
	nproc = args.nproc
	pool1 = multiprocessing.Pool ( processes=nproc )# len(num_list)) 
	pool_list = [pool1]
	rg = []
	
	for U in U_list:
		print (f"Diving into U = {U}...", flush=True)
		T_list = args.T

		# get num_list for each temperature 
		master_temperature_list     = []
		master_num_list      = []
		master_index_list    = []
		rg_dict    = {}
		ntraj_dict = {}
		for T in temperature_list:
			num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
			master_num_list.extend ( num_list )
			for num in num_list:
				master_index_list.append (get_starting_ind (U, T, num, dop, "energydump") )

			master_temperature_list.extend ( [T]*len( num_list ) )
			ntraj_dict[T] = len ( num_list )
			rg_dict[T]    = []

		# start multiprocessing... keeping in mind that each node only has 96 cores 
		# start splitting up master_num_list and master_temp_list 

		idx_range = len (master_num_list)//nproc + 1

		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_temperature_list[u_idx*nproc:], master_num_list [u_idx*nproc:], itertools.repeat (dop), itertools.repeat (coords_files), master_index_list[u_idx*nproc:] ) )
			else:
				results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_temperature_list[(u_idx)*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(coords_files), master_index_list[u_idx*nproc:(u_idx+1)*nproc] ) )

			print ("Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

			for k in range( len( master_temperature_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				rg_dict[master_temperature_list[u_idx*nproc + k]].append( results[k] )

		for temp in temperature_list:
			rg.extend( rg_dict[temp] )


	pool1.close()
	pool1.join()
	print (f"rg_max = {rg_max}...")
	i=0

	ax.hist(rg, density=True)
	# ax.set_xlim(1.0, 7.0)

	plt.savefig   ( args.pn, bbox_inches='tight', dpi=1200)
	##################################
	stop = time.time() 

	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

