#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np 
import re 
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import matplotlib.ticker as tck
import pandas as pd
import os
import time 
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/py_analysis')
import aux 
import multiprocessing 
import itertools
from sklearn.linear_model import LinearRegression 

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

'''
This code will take in a trajectory file generated by my MonteCarlo engine and 
gives you the flory exponent
'''
''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''


import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the flory exponent from that file.")
parser.add_argument('--integrated-database', dest='df', metavar='df', action='store', type=str, help='Name of dump file.')
parser.add_argument('--png-name', dest='pn', metavar='imagename', action='store', type=str, help='Name of image.')
args = parser.parse_args() 

# divnorm = matplotlib.colors.Normalize (vmin=-3, vmax=0)
grey_norm  = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.1, vmax=0 ) 
green_norm = matplotlib.colors.Normalize  ( vmin=-3.0, vmax=-0.1 ) 

def color_finder (H_mix):
	if H_mix <= -0.1:
		return cm.Greens_r (green_norm (H_mix) )
	elif H_mix >-0.1:
		return cm.Greys (grey_norm (H_mix) )

def chi_calculator (e_mm_a, e_mm_n, e_ms_a, e_ms_n, T):
	g = 0.25
	Zmm = g*np.exp (-1/T * e_mm_a) + (1-g)*np.exp(-1/T * e_mm_n)
	Zms = g*np.exp (-1/T * e_ms_a) + (1-g)*np.exp(-1/T * e_ms_n)
	fmm_a = g*np.exp (-1/T * e_mm_a) / Zmm 
	fms_a = g*np.exp (-1/T * e_ms_a) / Zms
	chi   = ( ( fms_a*e_ms_a + (1-fms_a)*e_ms_n ) - 0.5 * ( fmm_a*e_mm_a + (1-fmm_a)*e_mm_n) )/T
	
	return chi

def chi_extractor (U, T):
	E = aux.get_energy (U+"/geom_and_esurf.txt")
	chi = chi_calculator(E[0], E[1], E[2], E[3], T) 
	
	return chi

if __name__ == "__main__":
	start = time.time()
	##################################

	# temps = args.T # aux.dir2U ( os.listdir (".") )
	# U_list = aux.dir2U(os.listdir("."))
	fig = plt.figure   ( figsize=(4/1.6,3/1.6), constrained_layout=True )
	ax  = plt.axes() 
	plt.rcParams["axes.labelweight"] = "bold"
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	# ax.set (autoscale_on=False)
	# aux.gradient_image (ax, direction=0, extent=(0, 1, 0, 1), transform=ax.transAxes, cmap=plt.cm.coolwarm, cmap_range=(0.2, 0.8), alpha=1)
	i = 0

	##################################
	df = pd.read_csv (args.df, sep='|', header=0)
	U_list = np.unique(df["U"])
	print (U_list)
	i = 0
	chi_tot = []
	chi_list = []
	for U in U_list:
		chi_list.clear()
		# rgba_color = color_finder (U)
		nu = df.loc[df["U"] == U]
		for T in nu["T"]:
			chi_list.append ( chi_extractor (U, T) )
		chi_tot.extend(chi_list)
		# print (chi_list)
		# ax.errorbar(chi_list, nu["nu_mean"]/2, yerr=nu["nu_err"]/2, linewidth=0, capsize=2, \
		# ecolor='k', fmt='none', label='_nolegend_', marker='o', markersize=5)
		ax.plot(chi_list, nu["nu_mean"]/2, linewidth=0, marker='o',markersize=8/1.3, markeredgecolor='k', label="_nolabel_")# , c=rgba_color)
		i += 1
	df.insert (2, "chi", chi_tot, True)
	stop = time.time()
	ax.set_xscale ('symlog', linthresh=0.1)
	df.to_csv("INTEGRATED.csv", sep='|', float_format="%.2f", index=False)
	# yticks = np.arange(0.0, 0.9, 0.1) 
	# ax.set_yticks ( yticks )
	# ax.set_xticks (np.linspace(0, 1, 6))
	# ax.set_ylim   ( 0.0, 0.8 )
	# ax.set_xlim   ( -0.03, 1.03 )
	# ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	# ax.xaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	# ax.set_aspect('auto')
	# ax.set_xticklabels (ax.get_xticks(), weight='bold') 
	# ax.set_yticklabels (ax.get_yticks(), weight='bold') 
	# ax.yaxis.set_major_formatter(tck.StrMethodFormatter('{x:1.1f}') )
	# ax.xaxis.set_major_formatter(tck.StrMethodFormatter('{x:1.1f}') )
	plt.savefig   ( args.pn, dpi=1200)
	print ("Run time is {:.2f} seconds.".format(stop-start), flush=True)

