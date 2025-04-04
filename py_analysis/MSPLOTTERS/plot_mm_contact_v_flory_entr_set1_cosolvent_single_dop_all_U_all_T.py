#!/usr/licensed/anaconda3/2020.7/bin/python

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
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--path-to-excl', dest='pte', metavar='add/ress', action='store', type=str, help='Path to file.')
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = ["U10"] # aux.dir2U ( os.listdir(".") ) 
	plt.figure( figsize=(8,6) )
	
	PLOT_DICT = {}
	
	i=0
	ms_max = 208 # 25*2+(args.dop-2)*24
	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		ms_list = np.asarray([])
		ms_err  = np.asarray([])
		ms_mean = np.asarray([])
		temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
		# temperatures = [0.01, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0] # df[ df["U"] == U ]["T"]
		for temp in temperatures: 
			skip = 0
			ms_list = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
				ms_list = np.hstack ( (ms_list, np.mean(df["mm_aligned"].values[-2000:] ) ) )
                
			ms_err  = np.hstack ( (ms_err ,  (np.std (ms_list) / np.sqrt(len(num_list) ) ) ) )
			ms_mean = np.hstack ( (ms_mean,  np.mean (ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err )
		i += 1
		print("done!", flush=True)

	i=0
	chi_list = [0.1]# [0.1, 0.05, 0.01, 0.001, 0, -0.001, -0.01, -0.1, -0.2]
	for U in U_list:
		rgba_color = cm.PiYG(divnorm (chi_list[i]))
		df = pd.read_csv ("INTEGRATED-FLORY-EXPONENT-TYPE2.csv", sep='|')
		nu = df[df["U"]==U]
		nu = nu.loc[df["T"].isin(temperatures)]
		plt.errorbar ( nu["nu_mean"]/2, PLOT_DICT[U][0] / ms_max, yerr=PLOT_DICT[U][1]/ms_max, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     ( nu["nu_mean"]/2, PLOT_DICT[U][0] / ms_max, linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
		i += 1

	# plot excluded volume
	# contacts =  np.ones ( len(temperatures) ) 
	# if args.ev:
	# 	print ("In U = Uexcl...", flush=True)
	# 	df = pd.read_csv(args.pte+"/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], engine='python')
	# 	contacts = np.mean ( df["ms1_tot"].values + df["ms2_tot"].values ) * contacts
	# 	plt.errorbar ( temperatures, contacts/ms_max, yerr=0, fmt='-', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, markersize=10, linewidth=3)

	ax = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
	ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.set_xscale('log')
	ax.set_ylim((0.0 , 1.0))
	plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.2f}'))
	ax.set_yticks (np.linspace (0.00,1,6))
	ax.minorticks_on()
	ax.set_ylabel ("$C_{\\sigma}$")
	ax.set_xlabel ("Flory $\\nu$")
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)


