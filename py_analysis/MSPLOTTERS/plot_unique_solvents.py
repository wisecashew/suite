#!/home/satyend/.conda/envs/phase/bin/python

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
from pathlib import Path


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('--dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--U', dest='U', action='store', nargs='+', type=str, help='Name of forcefields.', default=[]) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--xllim',   dest='xllim', action='store', type=float, help="Lower x limit.", default=None)
parser.add_argument('--xulim',   dest='xulim', action='store', type=float, help="Upper x limit.", default=None)
parser.add_argument('--yllim',   dest='yllim', action='store', type=float, help="Lower y limit.", default=None)
parser.add_argument('--yulim',   dest='yulim', action='store', type=float, help="Upper y limit.", default=None)
parser.add_argument('--color', dest='col', metavar='color', nargs='+', action='store', type=str, help='Name of color of file.', default='coral') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
	# get the entire list of potential energy surfaces 
	fig = plt.figure(figsize=(1.65,1.65))
	ax  = plt.axes()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='y', labelsize=6)
	ax.tick_params(axis='x', labelsize=6)
	start = time.time()

	U_list = args.U
	
	PLOT_DICT = dict()
	
	i = 0
	Z = 26
	scale = args.dop * Z
	
	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		prop_list = np.asarray([])
		prop_err  = np.asarray([])
		prop_mean = np.asarray([])
		x_vals   = aux.dir2float(os.listdir(str(U)+"/DOP_"+str(args.dop)))
		for x in x_vals:
			skip = 0
			prop_list = np.asarray ([])
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(frac) ) ) )

			for num in num_list:
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(x)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "timestep"], engine='python', skiprows=skip)
				prop_list = np.hstack ((ms_list, np.mean(df["total"].values[-args.s:])))

			prop_err  = np.hstack ( (prop_err ,  (np.std (prop_list) / np.sqrt(len(num_list) ) ) ) )
			prop_mean = np.hstack ( (prop_mean,  np.mean (prop_list) ) )

		PLOT_DICT[U] = ( prop_mean, prop_err )
		i += 1
		print("done!", flush=True)

	i=0
	for U in U_list:
		plt.errorbar (np.array(x_vals), PLOT_DICT[U][0]/scale, yerr=PLOT_DICT[U][1]/scale, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     (np.array(x_vals), PLOT_DICT[U][0]/scale, linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, markersize=8/1.3, zorder=10, clip_on=False, c=args.col[i])
		i += 1

	ax.set_xlim(args.xllim, args.xulim)
	ax.set_ylim(args.yllim, args.yulim)
	ax.set_xticklabels([])
	ax.set_yticklabels([])  
	ax.minorticks_on()
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)


