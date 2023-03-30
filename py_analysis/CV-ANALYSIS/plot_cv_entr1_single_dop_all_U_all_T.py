#!/usr/licensed/anaconda3/2020.7/bin/python

import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/py_analysis")

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = ["U10"] # aux.dir2U ( os.listdir(".") )
	plt.figure( figsize=(8,6) )

	PLOT_DICT = {}

	i=0
	Tmax = []
	ms_max = 25*2+(args.dop-2)*24

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		ms_list = np.asarray([])
		ms_err  = np.asarray([])
		ms_mean = np.asarray([])
		temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )

		for temp in temperatures: 
			skip = 0
			ms_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
				f = df["energy"].values[-2000:]
				ms_list = np.hstack ( (ms_list, ( np.mean( f**2 ) - np.mean (f)**2 ) ) ) # / temp**2  ) )
			ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(len(num_list)) ) ) )
			ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err ) 
		i += 1
		# mm_max.append( np.max(ms_list) ) 
		print("done!", flush=True)

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
		chi_1 = chi_list[i] # aux.get_chi_cosolvent ( str(U)+"/geom_and_esurf.txt" )[args.cs]
		print ("chi_1 = ", chi_1)
		rgba_color = cm.PiYG(divnorm (chi_1))
		plt.errorbar ( temperatures, PLOT_DICT[U][0]/(args.dop), yerr=PLOT_DICT[U][1]/(args.dop), linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     ( temperatures, PLOT_DICT[U][0]/(args.dop), linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)

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
	contacts =  np.ones ( len(temperatures) ) 

	ax = plt.axes ()
	ax.tick_params ( direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
	ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )

	ax.set_xscale('log')
	plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.1f}'))
	ax.minorticks_on()
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)

