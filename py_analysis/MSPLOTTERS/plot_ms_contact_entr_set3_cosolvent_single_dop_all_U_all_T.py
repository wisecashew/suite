#!/usr/licensed/anaconda3/2020.7/bin/python

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


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(4/1.6,3/1.6), constrained_layout=True )
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=9, pad=5)
	ax.tick_params(axis='y', labelsize=9)
	norm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 ) # this is for entropy 
	font = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size':11}

	color1 = np.array([131, 159, 192]) / 255.0  # #839FC0 in RGB
	color2 = np.array([241, 156, 118]) / 255.0   # #ED8151 in RGB
	cmap = colors.LinearSegmentedColormap.from_list('custom', [color1, "white", color2])
	cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)
	y = np.array([0,1])

	df = pd.read_csv ("INTEGRATED-FLORY-EXPONENT-TYPE2.csv", sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
	df = df.loc[df["U"] == "U4"]
	temperatures = df["T"].values
	# df = df[df["T"].isin(temperatures)]
	print (df)

	x_old = temperatures
	y_old = df["nu_mean"].values/2

	x_pred = np.logspace (-2, 2, 10000)
	y_pred = np.interp (x_pred, x_old, y_old)

	col_dict = dict()

	for i in range (len(x_pred)):
		col_dict[ x_pred[i] ] = y_pred[i]

	X, Y = np.meshgrid (x_pred, y)
	Z = np.zeros ((2, len(x_pred)))

	for i in range (len(x_pred)):
		for j in range (2):
			Z[j, i]= col_dict [x_pred[i]]

	ax.pcolormesh (X, Y, Z, cmap=cmap, norm=cnorm, shading="auto")

	start = time.time()

	U_list = ["U4"] # aux.dir2U ( os.listdir(".") ) 
	
	PLOT_DICT = dict()
	
	i=0
	ms_max = 25*2+(args.dop-2)*24
	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		ms_list = np.asarray([])
		ms_err  = np.asarray([])
		ms_mean = np.asarray([])
		temperatures = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
		for temp in temperatures:
			skip = 0
			ms_list = np.asarray ([])
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
				ms_list = np.hstack ( (ms_list, np.mean(df["ms1_aligned"].values[-2000:] ) ) )

			ms_err  = np.hstack ( (ms_err ,  (np.std (ms_list) / np.sqrt(len(num_list) ) ) ) )
			ms_mean = np.hstack ( (ms_mean,  np.mean (ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err )
		i += 1
		print("done!", flush=True)

	i=0
	chi_list = [0.1]
	for U in U_list:
		rgba_color = cm.PiYG(divnorm (chi_list[i]))
		idx = [0, 1, 3, 5, 6, 7, 9, 10, 11, 13, 15]
		plt.errorbar ( np.array(temperatures), PLOT_DICT[U][0] / ms_max, yerr=PLOT_DICT[U][1]/ms_max, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     ( np.array(temperatures), PLOT_DICT[U][0] / ms_max, linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10, zorder=10, clip_on=False)
		i += 1

	ax.set_xscale('log')
	yticks = np.arange(0.0, 1.2, 0.2)
	ax.set_yticks ( yticks )
	ax.set_yticklabels (ax.get_yticks(), fontdict=font) 
	ax.set_ylim   ( 0.0, 1.0 )
	ax.set_xlim   ( 0.01, 100 )
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticklabels (["$10^{-2}$", "$10^{-1}$", "$10^0$", "$10^1$", "$10^2$"], fontdict=font)
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
	ax.yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}') )
	ax.set_aspect ('auto')
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)


