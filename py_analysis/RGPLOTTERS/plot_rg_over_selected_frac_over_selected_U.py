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
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-dop',           metavar='DOP',   dest='dop',            type=int,      action='store', help='enter a degree of polymerization.')
parser.add_argument('-s',             metavar='S',     type=int,              dest='s',      action='store', help='start parsing after this move number (not index or line number in file).', default=100)
parser.add_argument('--yllim',        metavar='yllim', type=float,           dest='yllim',   action='store', help='enter lower ylimit.', default=0)
parser.add_argument('--yulim',        metavar='yulim', type=float,           dest='yulim',   action='store', help='enter upper ylimit.', default=10)
parser.add_argument('--frac',         metavar='frac', dest='frac',           action='store', nargs='+',      type=float, help='Enter fractions you want probed.')
parser.add_argument('--U',            metavar='U',    dest='U',              action='store', nargs='+',      type=str,   help='Enter U you want probed.')
parser.add_argument('-nproc',         metavar='N',    type=int,              dest='nproc',   action='store', help='Request these many proccesses.')
parser.add_argument('--coords',       dest='c',       metavar='coords.txt',  action='store', type=str,       help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--show-legends', dest='sl',      action='store_true',   help='Name of energy dump file to parse information.', default=False)
parser.add_argument('--show-yticks',  dest='syticks', action='store_true',   help='Show y ticks on plot.', default=False)
parser.add_argument('--png-name',     dest='pn',      metavar='imagename',   action='store', type=str,       help='Name of image file', default='rg_plot')
args = parser.parse_args() 

def get_starting_ind ( U, frac, num, dop, dumpfile):
    filename = U + "/DOP_" + str(dop) + "/" + str(frac) + "/" + dumpfile + "_" + str(num) + ".mc"
    df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    L = len(df["energy"])
    return int(df["time_step"].values[L-2000])


def infiltrate_coords_get_rg ( U, frac, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(frac) + "/"+ coords_files + "_" + str(num)+".mc" 
    edge = aux.edge_length (dop)
    master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
    rg = aux.get_Rg(master_dict, edge, edge, edge) 
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
	cols    = color_map("coral", "darkred", len(U_list))
	print (U_list, flush=True)
	PLOT_DICT = {} 
	dop            = args.dop
	coords_files   = args.c
	starting_index = args.s

######
	fig = plt.figure( figsize=(2.5,2.5), constrained_layout=True )
	ax  = plt.axes() 
# plt.rcParams["axes.labelweight"] = "bold"
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	i = 0 

	rg_max = 1 # (1+np.sqrt(2)+np.sqrt(3))/(3*6**0.5) * (dop**0.57) 
# instantiating pool
	nproc = args.nproc
	pool1 = multiprocessing.Pool ( processes=nproc )# len(num_list)) 

	pool_list = [pool1] # , pool2]

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
			num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(f) ) ) ) )
			master_num_list.extend ( num_list )
			for num in num_list:
				master_index_list.append (get_starting_ind (U, f, num, dop, "energydump") )

			master_frac_list.extend ( [f]*len( num_list ) )
			ntraj_dict[f] = len ( num_list )
			rg_dict[f]    = []

        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 

		idx_range = len (master_num_list)//nproc + 1

		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_frac_list[u_idx*nproc:], master_num_list [u_idx*nproc:], itertools.repeat (dop), itertools.repeat (coords_files), master_index_list[u_idx*nproc:] ) )
			else:
				results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_frac_list[(u_idx)*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(coords_files), master_index_list[u_idx*nproc:(u_idx+1)*nproc] ) )

			print ("\tPool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

			for k in range( len( master_frac_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				rg_dict[master_frac_list[u_idx*nproc + k]].append( results[k] )

		for f in frac_list:
			rg_mean.append( np.mean ( rg_dict[f] ) )
			rg_std.append ( np.std  ( rg_dict[f] ) / np.sqrt(ntraj_dict[f] ) )

		PLOT_DICT [U] = ( np.asarray(rg_mean), np.asarray(rg_std) )


	pool1.close()
	pool1.join()
	print ("rg_max = ", rg_max)
	i=0

	for idx,U in enumerate(U_list):
		ax.errorbar ( frac_list, PLOT_DICT[U][0]/rg_max, yerr= PLOT_DICT[U][1]/rg_max, linewidth=1, capsize=2, color='k', fmt='none', label='_nolegend_')
		ax.plot     ( frac_list, PLOT_DICT[U][0]/rg_max, marker='o', markeredgecolor='k', linestyle='-', c=cols[idx], linewidth=1, markersize=8/1.3, label=f'{U}', clip_on=False, zorder=10) 
		i += 1

	#########
	# plt.rcParams['text.usetex'] = True
	ax.minorticks_on()
	ax.set_xlim((0.0, 1.0))
	ax.set_ylim((args.yllim, args.yulim))
	ax.set_xticks (np.linspace (0, 1, 6))
	ax.set_xticklabels ([])
	if args.syticks:
		plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
	else:
		ax.set_yticklabels ([])

	ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator())

	if args.sl:
		ax.legend(loc="upper right", fontsize=6)
	else:
		print ("No legends in image.")
		ax.legend().remove()
	plt.savefig   ( args.pn, bbox_inches='tight', dpi=1200)

##################################
	stop = time.time() 

	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

