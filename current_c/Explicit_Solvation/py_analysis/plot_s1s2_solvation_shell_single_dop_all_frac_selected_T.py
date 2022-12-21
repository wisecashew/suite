#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux 
import os 
import re
import multiprocessing
import itertools

parser = argparse.ArgumentParser(description="Read a lattice dump file and obtain the aligned s1s2 contacts in the solvation shell.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-T', metavar='T', dest='T', action='store', nargs='+', type=float, help='Enter temperatures you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('-ymax', metavar='ymax', type=float, dest='ymax', action='store', default=2.8, help='Request upper limit of plot.')
parser.add_argument('--png-name', dest='pn', metavar='imagename', action='store', type=str, help='Name of image file', default='rg_plot')
args = parser.parse_args() 
# divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 
def mysorter_f (x):
	new_str = ""
	for m in x:
		if m.isdigit():
			new_str = new_str+m
	return float(new_str)

if __name__=="__main__":

	U_list = aux.dir2U (os.listdir("."))
	PLOT_DICT = {}
	dop = args.dop
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort() 

	fig = plt.figure ( figsize=(8,6) )
	ax  = plt.axes () 
	
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='x', labelsize=16)
	ax.tick_params(axis='y', labelsize=16)
	i = 0 
	nproc = args.nproc
	pool1 = multiprocessing.Pool ( processes=nproc )# len(num_list)) 


	for T in temperatures:
		frac_list = [] 
		print("Inside T = " + str(T) + ", and N = " + str(dop) + "...", flush=True )
		aligned_ave = [] 
		aligned_err = []

		master_U_list = []
		master_num_list = []
		aligned_dict = {}
		ntraj_dict   = {}

		for U in U_list:
			frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt") )
			num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) ) 
			master_num_list.extend ( num_list ) 
			master_U_list.extend ( [U]*len(num_list) ) 
			ntraj_dict[U] = len(num_list)
			aligned_dict[U] = []

		idx_range = len (master_num_list)//nproc + 1
		print ("Number of loops: {}".format (idx_range), flush=True)
		for u_idx in range (idx_range):
			print ("Loop index: {}".format (u_idx), flush=True)
			if u_idx == idx_range-1:
				results = pool1.starmap ( aux.extract_s1s2_aligned_solvation_shell, zip( master_U_list[u_idx*nproc:], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:] ) )
			else:
				results = pool1.starmap ( aux.extract_s1s2_aligned_solvation_shell, zip( master_U_list[u_idx*(nproc):(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc] ) )

		print ("Poolhas been closed. This pool had {} threads.".format(len(results)), flush=True)
		for k in range ( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
			aligned_dict[master_U_list[u_idx*nproc + k]].append ( results[k] )

		sorted_U_list = list (np.unique (master_U_list)) 
		sorted_U_list.sort (key=mysorter_f)

		for U in sorted_U_list:
			aligned_ave.append ( np.mean (aligned_dict[U]) )
			aligned_err.append ( np.std  (aligned_dict[U]) / np.sqrt (master_U_list.count(U) ) )
		
		print (aligned_ave)
		PLOT_DICT[T] = (np.array(aligned_ave), np.array(aligned_err) )


	pool1.close()
	pool1.join()


	i=0
	if len(temperatures) == 1:
		colors = cm.coolwarm (np.linspace(0,1,3)) 
		i = 1
	else:
		colors = cm.coolwarm (np.linspace(0,1,len(temperatures)))
	for T in temperatures:
		ax.errorbar ( frac_list, PLOT_DICT[T][0], yerr= PLOT_DICT[T][1], linewidth=1, capsize=2, color='k', fmt='none', label='_nolegend_')
		ax.plot     ( frac_list, PLOT_DICT[T][0], marker='o', markeredgecolor='k', linestyle='-', linewidth=3, c=colors[i], label='_nolegend_', markersize=10 ) 
		i += 1



	my_cmap = cm.coolwarm
	sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0) ) 
	ax = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )

	cbar = plt.colorbar(sm, orientation='vertical')
	cbar.set_ticks( [] )
	# cbar.ax.tick_params (labelsize=14)
	cbar.set_ticklabels( [] )
	# cbar.ax.set_ylabel ( "Strength of better solvent \n", fontsize=16, rotation=270 ) 
	# ax.set_xscale('log')
	# ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
	# ax.yaxis.get_major_locator().set_params(integer=True)
	ax.set_ylim(bottom=0.5)
	ax.set_yticks (np.linspace (0.6, 2.4, 10))
	#yticks = np.linspace(0.2, 0.8, 7) 
	#yticks = np.hstack ((yticks, 0.57))
	#yticks[0] = 0.25 
	#yticks[-2] = 0.75
	#yticks = np.hstack((yticks, 0.33))
	ax.set_xticks (np.linspace (0, 1, 6) )
	ax.minorticks_on()
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	#ax.set_yticks ( yticks ) 
	
	# plt.legend(U_list)
	# plt.savefig   (args.pn + ".png", dpi=1000)


