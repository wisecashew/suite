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
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--U', dest='U', action='store', nargs='+', type=str, help='Name of forcefields.', default=[]) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(2,2), constrained_layout=True )
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='y', labelsize=6)
	ax.tick_params(axis='x', labelsize=6)
	start = time.time()

	U_list = args.U
	
	PLOT_DICT = dict()
	
	i=0
	ms_max = 208 # 25*2+(args.dop-2)*24
	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		ms_list = np.asarray([])
		ms_err  = np.asarray([])
		ms_mean = np.asarray([])
		fracs   = aux.dir2float(os.listdir(str(U)+"/DOP_32"))
		# temperatures = [0.01, 0.02, 0.03, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.9, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
		for f in fracs:
			skip = 0
			ms_list = np.asarray ([])
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(f) ) ) )

			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(f)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
				ms_list = np.hstack ( (ms_list, np.mean(df["mm_tot"].values[-args.s:] ) ) )

			ms_err  = np.hstack ( (ms_err ,  (np.std (ms_list) / np.sqrt(len(num_list) ) ) ) )
			ms_mean = np.hstack ( (ms_mean,  np.mean (ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err )
		i += 1
		print("done!", flush=True)

	i=0
	chi_list = [-0.2]
	for U in U_list:
		rgba_color = "#B91F72"
		plt.errorbar ( np.array(fracs), PLOT_DICT[U][0], yerr=PLOT_DICT[U][1], linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     ( np.array(fracs), PLOT_DICT[U][0], linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, label=f"U = {U}", markersize=8, zorder=10, clip_on=False)
		i += 1

	ax.legend(loc="upper right", prop={'size':5})
	# yticks = np.arange(0.0, 1.2, 0.2)
	# ax.set_yticks ( yticks )
	# ax.set_yticklabels ([]) # ax.get_yticks(), fontdict=fdict, font=fpath)
	ax.set_ylim   ( 0.0, 210 )
	# ax.set_xlim   ( 0.0, 1.0 )
	# ax.set_xticks ([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.set_xticklabels ([0.01, 0.1, 1.0, 10.0, 100.0], fontdict=fdict, font=fpath)
	ax.set_xticks ( np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	# ax.set_xticklabels ([])
	# ax.yaxis.set_minor_locator (tck.AutoMinorLocator())
	# ax.yaxis.set_major_formatter(tck.StrMethodFormatter('{x:1.1f}') )
	# ax.set_yticklabels ([])
	ax.set_aspect ('auto')
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)


