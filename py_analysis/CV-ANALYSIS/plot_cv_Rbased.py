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

parser = argparse.ArgumentParser(description="Get the heat capacities from our simulation.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.', default=32) 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--yulim', dest='yulim', action='store', type=float, help='Provide an upper limit to y-axis.', default=None)
parser.add_argument('--yllim', dest='yllim', action='store', type=float, help='Provide a lower limit to y-axis.', default=None)
parser.add_argument('--T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature list.')
parser.add_argument('--database', dest='db', action='store', type=str, help='The database where the flory exponents is stored.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump')
parser.add_argument('--fluc-style', dest='fluc_style', action='store', type=str, help='Type of fluctuations to measure: total, monomer-monomer, monomer-solvent.')
parser.add_argument('--ylabels', dest='ylabels', action='store_true', default=False, help='Provide this label if you want to see ylabels.')
parser.add_argument('--no-kt2', dest='nokt2', action='store_true', default=False, help='Provide this label if you want to scale by kT^2.')
parser.add_argument('--R', dest='R', metavar='RX', action='store', type=str, help='Name of forcefield in database to plot.')
parser.add_argument('--U', dest='U', metavar='UX', action='store', type=str, help='Name of forcefield directory.')
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()
nokt2 = args.nokt2
divnorm    = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy
color_dict = {"R0":'#369DE8', "R1": '#1FB967', "R2":'#B9B41F', "R3":'#B91F72'}

def total_heat_capacity (U_list, PLOT_DICT):

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for temp in temperatures: 
			skip = 0
			energy_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = df["energy"].values[-2000:]/args.dop
				if nokt2:
					energy_list = np.hstack ( (energy_list, np.mean( f**2 ) - np.mean (f)**2) )
				else:
					energy_list = np.hstack ( (energy_list, ( np.mean( f**2 ) - np.mean (f)**2 ) / temp**2 ) )
			energy_err  = np.hstack ( (energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)) ) ) )
			energy_mean = np.hstack ( (energy_mean, np.mean(energy_list) ) )

		PLOT_DICT[U] = ( energy_mean, energy_err )
		print (f"done with computations for {U}!", flush=True)
	
	return

def mm_heat_capacity (U_list, PLOT_DICT):

	for U in U_list:
		path_to_top = U + "/geom_and_esurf.txt"
		energetic_params = aux.get_energy_target (path_to_top)
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for temp in temperatures: 
			skip = 0
			energy_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = (df["mm_aligned"].values[-2000:] * energetic_params[0] + df["mm_naligned"].values[-2000:] * energetic_params[1])/args.dop
				if nokt2:
					energy_list = np.hstack ( (energy_list, np.mean( f**2 ) - np.mean (f)**2) )
				else:
					energy_list = np.hstack ( (energy_list, ( np.mean( f**2 ) - np.mean (f)**2 ) / temp**2 ) )
			energy_err  = np.hstack ( (energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)) ) ) )
			energy_mean = np.hstack ( (energy_mean, np.mean(energy_list) ) )

		PLOT_DICT[U] = ( energy_mean, energy_err )
		print (f"done with computations for {U}!", flush=True)
	
	return

def ms_heat_capacity (U_list, PLOT_DICT):

	for U in U_list:
		path_to_top = U + "/geom_and_esurf.txt"
		energetic_params = aux.get_energy_target (path_to_top)
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for temp in temperatures: 
			skip = 0
			energy_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = (df["ms1_aligned"].values[-2000:] * energetic_params[2] + df["ms1_naligned"].values[-2000:] * energetic_params[3])/args.dop
				if nokt2:
					energy_list = np.hstack ( (energy_list, np.mean( f**2 ) - np.mean (f)**2) )
				else:
					energy_list = np.hstack ( (energy_list, ( np.mean( f**2 ) - np.mean (f)**2 ) / temp**2 ) )
			energy_err  = np.hstack ( (energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)) ) ) )
			energy_mean = np.hstack ( (energy_mean, np.mean(energy_list) ) )

		PLOT_DICT[U] = ( energy_mean, energy_err )
		print (f"done with computations for {U}!", flush=True)
	
	return


def cross_heat_capacity (U_list, PLOT_DICT):

	for U in U_list:
		path_to_top = U + "/geom_and_esurf.txt"
		energetic_params = aux.get_energy_target (path_to_top)
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for temp in temperatures: 
			skip = 0
			energy_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f_ms = (df["ms1_aligned"].values[-2000:] * energetic_params[2] + df["ms1_naligned"].values[-2000:] * energetic_params[3])/args.dop
				f_mm = (df["mm_aligned"].values[-2000:] * energetic_params[0] + df["mm_naligned"].values[-2000:] * energetic_params[1])/args.dop
				if nokt2:
					energy_list = np.hstack ((energy_list, np.mean(f_mm*f_ms) - np.mean(f_mm)*np.mean(f_ms)))
				else:
					energy_list = np.hstack ((energy_list, (np.mean(f_mm*f_ms) - np.mean(f_mm)*np.mean(f_ms)) / temp**2))
			energy_err  = np.hstack ( (energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)) ) ) )
			energy_mean = np.hstack ( (energy_mean, np.mean(energy_list) ) )

		PLOT_DICT[U] = ( energy_mean, energy_err )
		print (f"done with computations for {U}!", flush=True)
	
	return



if __name__=="__main__":

	temperatures = args.T
	print (f"Temperatures to be probed = {args.T}...")

	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(1.7,1.7) )
	lsize = 14
	fig.tight_layout()
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	norm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 )
	fdict = {'color': 'black', 'weight': 'normal', 'size': lsize}

	# plotted the background

	PLOT_DICT = {}

	i=0

	U_list = [args.U]

	fluc_style = args.fluc_style

	if fluc_style == "total":
		total_heat_capacity (U_list, PLOT_DICT)

	elif fluc_style == "monomer-monomer":
		mm_heat_capacity (U_list, PLOT_DICT)

	elif fluc_style == "monomer-solvent":
		ms_heat_capacity (U_list, PLOT_DICT)

	elif fluc_style == "cross":
		cross_heat_capacity (U_list, PLOT_DICT)

	else:
		print(f"There is no fluctuation with the label \'{fluc_style}\'. Exiting...")
		exit()

	i=0

	for U in U_list:
		rgba_color = color_dict[args.R]
		ax.errorbar ( temperatures, PLOT_DICT[U][0], yerr=PLOT_DICT[U][1], linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		ax.plot     ( temperatures, PLOT_DICT[U][0], linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, color=rgba_color, label="_nolabel_", markersize=8/1.3, clip_on=False, zorder=10)

	# plot the background
	color1 = np.array([131, 159, 192]) / 255.0  # #839FC0 in RGB
	color2 = np.array([137, 245, 162]) / 255.0   # #ED8151 in RGB
	cmap = colors.LinearSegmentedColormap.from_list('custom', [color2, "white", color1])
	cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)

	if args.yulim == None:
		yulim = np.max(PLOT_DICT[U][0]) * 1.1 if np.max(PLOT_DICT[U][0]) > 0 else np.max(PLOT_DICT[U][0])*0.9
	else:
		yulim = args.yulim

	if args.yllim == None:
		yllim = np.min(PLOT_DICT[U][0]) * 0.9 if np.max(PLOT_DICT[U][0]) > 0 else np.min(PLOT_DICT[U][0]) * 1.1
	else:
		yllim = args.yllim

	y = np.array([yllim, yulim])

	df = pd.read_csv (args.db, sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
	df = df.loc[df["U"] == args.R]
	temperatures = df["T"].values

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

	if args.ylabels:
		pass
	else:
		ax.set_yticklabels ([])

	ax.set_ylim   (yllim, yulim)
	ax.set_yticks (np.linspace(yllim, yulim, 5))
	ax.set_xlim   (0.01, 100)
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticklabels ([])
	ax.set_xscale ("log")
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticks (np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	ax.set_xticklabels ([])
	ax.set_aspect ('auto')
	plt.savefig   (args.pn, bbox_inches='tight', dpi=1200)

