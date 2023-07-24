#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np 
from pathlib import Path
import re 
import matplotlib
from matplotlib import rcParams
import matplotlib.style as mpl
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.colors as colors
import pandas as pd
import os
import time 
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux 
import multiprocessing 
import itertools
from sklearn.linear_model import LinearRegression 
import argparse

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = ["U9"] # aux.dir2U ( os.listdir(".") )
	fig = plt.figure( figsize=(1.7,1.7), constrained_layout=True )
	ax  = plt.axes  ()

	ax.tick_params ( direction='in', bottom=True, top=True, left=True, right=True, which='both')
	
	# Define the gradient colors
	color1 = np.array([131, 159, 192]) / 255.0  # #839FC0 in RGB
	color2 = np.array([137, 245, 162]) / 255.0   # #ED8151 in RGB
	cmap = colors.LinearSegmentedColormap.from_list('custom', [color2, "white", color1])
	cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)
	y = np.array([0,0.01])

	df = pd.read_csv ("INTEGRATED-FLORY-EXPONENT-TYPE2.csv", sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
	df = df.loc[df["U"] == "U9"]
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
	

	PLOT_DICT = {}

	i=0
	Tmax = []
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
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = df["energy"].values[-2000:]/args.dop
				ms_list = np.hstack ( (ms_list, ( np.mean( f**2 ) - np.mean (f)**2 )/temp**2 ) )
			ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(len(num_list)) ) ) )
			ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err )
		i += 1
		print ("done!", flush=True)

	# get the maximum value of fluctuation
	for key in PLOT_DICT:
		if ms_max < np.max(PLOT_DICT[key][0]):
			ms_max = np.max(PLOT_DICT[key][0])

	print ("ms_max = ",ms_max)
	# ms_max = 1
	i=0
	ymin = 1
	chi_list = [0.1, 0.05, 0.01, 0.001, 0, -0.001, -0.01, -0.1, -0.2]
	f = open ("CV-output.mc", 'w')
	for U in U_list:
		chi_1 = chi_list[i]
		print ("chi_1 = ", chi_1)
		rgba_color = "#B91F72" # cm.PiYG(divnorm (chi_list[i]))
		ax.errorbar ( temperatures, PLOT_DICT[U][0]/(args.dop), yerr=PLOT_DICT[U][1]/(args.dop), linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		ax.plot     ( temperatures, PLOT_DICT[U][0]/(args.dop), linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, color=rgba_color, label="_nolabel_", markersize=8/1.3, clip_on=False, zorder=10)

		i += 1
		f.write ("U = "+U+":\n")
		f.write ("T:  ")
		for t in temperatures:
			f.write ("{:>2.2f} ".format(t))
			f.write("\n")
		f.write ("Cv: ")
		h = 0
		for c in PLOT_DICT[U][0]:
			f.write ("{:>2.2e} ".format(c/temperatures[h]**2)) 
			h += 1
		f.write ("\n")

	f.close()
	# plot excluded volume
	yticks = np.arange(0.0, 0.012, 0.002)
	ax.set_yticks ( yticks )
	ax.set_ylim   ( 0.0, 0.01 )
	ax.set_xlim   ( 0.01, 100 )
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticklabels ([]) # [0.01, 0.1, 1.0, 10.0, 100.0], fontdict=fdict, font=fpath)
	ax.set_xscale ("log")
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticks ( np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	# ax.yaxis.set_major_formatter(tck.StrMethodFormatter('{x:1.1f}') )
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	# ax.set_yticklabels ([])
	ax.set_xticklabels ([])
	ax.set_aspect ('auto')

	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)

