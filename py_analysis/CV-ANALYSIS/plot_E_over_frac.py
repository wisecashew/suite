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
parser.add_argument('-dop'     , dest='dop'    , action='store', type=int, help='Provide degree of polymerization.', default=32) 
parser.add_argument('-s'       , dest='s'      , action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--figsize', dest='figsize', action='store', nargs=2 , type=float, help='Enter dimensions of image (width, height) (default: (2.5, 2.5)', default=[2.5,2.5])
parser.add_argument('--nyticks'  , dest='nyticks', action='store', type=int, help='Provide number of major yticks.', default=5)
parser.add_argument('--yulim'    , dest='yulim', action='store', type=float, help='Provide an upper limit to y-axis.', default=None)
parser.add_argument('--yllim'    , dest='yllim', action='store', type=float, help='Provide a lower limit to y-axis.', default=None)
parser.add_argument('--temp'     , dest='temp', action='store', type=float, help='Provide a temperature.',default=None)
parser.add_argument('--frac'     , dest='frac', action='store', type=float, nargs='+', help='Provide a frac list.',default=None)
parser.add_argument('--database' , dest='db', action='store', type=str, help='The database where the flory exponents is stored.', default=None)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump')
parser.add_argument('--fluc-style', dest='fluc_style', action='store', type=str, help='Type of fluctuations to measure: total, chain.')
parser.add_argument('--ylabels', dest='ylabels', action='store_true', default=False, help='Provide this label if you want to see ylabels.')
parser.add_argument('--no-kt2', dest='nokt2', action='store_true', default=False, help='Provide this label if you want to scale by kT^2.')
parser.add_argument('--R', dest='R', metavar='RX', action='store', type=str, help='Name of forcefield in database to plot.', default=None)
parser.add_argument('--color', dest='cols', metavar='color', action='store', nargs='+', type=str, help='Enter color of plot.', default=None)
parser.add_argument('--U', dest='U', metavar='UX', action='store', nargs='+', type=str, help='Name of forcefield directory.')
parser.add_argument('--clip-on', dest='clipon', action='store_true', help='Option to clip points mpl.', default=False)
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()
nokt2 = args.nokt2
divnorm    = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy
color_dict = {"R0":'#369DE8', "R1": '#1FB967', "R2":'#B9B41F', "R3":'#B91F72'}

def total_energy (U_list, PLOT_DICT, frac_list):

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for frac in frac_list:
			skip = 0
			energy_list = np.asarray([])
			num_list    = np.unique(aux.dir2nsim(os.listdir(str(U)+"/DOP_"+str(args.dop)+"/"+str(frac))))
			for num in num_list:
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(frac)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = df["energy"].values[-2000:]/args.dop
				energy_list = np.hstack((energy_list, np.mean(f)))
			energy_err  = np.hstack ((energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)))))
			energy_mean = np.hstack ((energy_mean, np.mean(energy_list)))

		PLOT_DICT[U] = (energy_mean, energy_err)
		print (f"done with computations for {U}!", flush=True)

	return


def total_chain_energy (U_list, PLOT_DICT, frac_list):

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		energy_list = np.asarray([])
		energy_err  = np.asarray([])
		energy_mean = np.asarray([])

		for frac in frac_list:
			skip = 0
			energy_list = np.asarray([])
			num_list    = np.unique(aux.dir2nsim(os.listdir(str(U)+"/DOP_"+str(args.dop)+"/"+str(frac))))
			energy      = aux.get_energy_target(str(U)+"/DOP_"+str(args.dop)+"/"+str(frac)+"/"+"geom_and_esurf.txt")
			for num in num_list:
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(frac)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = energy[0]*df["mm_aligned"].values + energy[1]*df["mm_naligned"].values + energy[2]*df["ms1_aligned"].values + energy[3]*df["ms1_naligned"].values + energy[4]*df["ms2_aligned"].values + energy[5]*df["ms2_naligned"].values
				energy_list = np.hstack((energy_list, np.mean(f)/args.dop))
			energy_err  = np.hstack ((energy_err,  (np.std(energy_list) / np.sqrt(len(num_list)))))
			energy_mean = np.hstack ((energy_mean, np.mean(energy_list)))

		PLOT_DICT[U] = (energy_mean, energy_err)
		print (f"done with computations for {U}!", flush=True)

	return



if __name__=="__main__":

	frac = args.frac
	temp = args.temp
	print (f"Fraction to be probed = {frac}...")
	fsize = args.figsize

	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(fsize[0],fsize[1]) )
	lsize = 14
	fig.tight_layout()
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	norm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 )
	fdict = {'color': 'black', 'weight': 'normal', 'size': lsize}

	# plotted the background

	PLOT_DICT = {}

	i=0

	U_list = args.U

	fluc_style = args.fluc_style

	if fluc_style == "total":
		total_energy (U_list, PLOT_DICT, frac)

	elif fluc_style == "chain":
		total_chain_energy (U_list, PLOT_DICT, frac)

	else:
		print(f"There is no fluctuation with the label \'{fluc_style}\'. Exiting...")
		exit()

	for idx, U in enumerate(U_list):
		ax.errorbar (frac, PLOT_DICT[U][0], yerr=PLOT_DICT[U][1], linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		ax.plot     (frac, PLOT_DICT[U][0], linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, color=args.cols[idx], label="_nolabel_", markersize=8/1.3, clip_on=args.clipon, zorder=10)

	if args.ylabels:
		pass
	else:
		ax.set_yticklabels ([])
	yllim = args.yllim
	yulim = args.yulim

	if yllim is None or yulim is None:
		pass
	else:
		ax.set_ylim   (yllim, yulim)
		ax.set_yticks (np.linspace(yllim, yulim, args.nyticks))
	ax.set_xlim   (0, 1)
	ax.set_xticks (np.linspace(0, 1, 6))
	ax.set_xticklabels ([])
	# ax.set_xticks (np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	ax.minorticks_on()
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	ax.set_xticklabels ([])
	ax.set_aspect ('auto')
	plt.savefig   (args.pn, bbox_inches='tight', dpi=1200)

