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

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
# parser.add_argument('--dump-file', dest='e', metavar='FLORY_DATA.mc', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
# parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
# parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.')

args = parser.parse_args()

# divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

if __name__=="__main__":

	results = aux.extract_s1s2_aligned_solvation_shell ("U11", 32, 0.01, 1)

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


