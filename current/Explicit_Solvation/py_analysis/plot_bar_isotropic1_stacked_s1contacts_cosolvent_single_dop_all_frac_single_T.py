#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
from matplotlib.ticker import StrMethodFormatter
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as tck
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
# parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
# parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.')
# parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
# parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(8,6))
	file_list = ["s1s1.dat", "s1s2.dat", "m1sn.dat"] 
	for f in file_list:
		if f == "s1s1.dat":
			df1 = pd.read_csv(f, sep='\t', engine='python', names=["frac", "contacts"], skiprows=1)
			ax.bar(df1["frac"].values, df1["contacts"].values, color='crimson', width=0.05, edgecolor='k') 
		elif f == "s1s2.dat":
			df2 = pd.read_csv(f, sep='\t', engine='python', names=["frac", "s1s2", "s2s1"], skiprows=1)
			print (df2)
			ax.bar(df2["frac"].values, df2["s1s2"].values, color='goldenrod', width=0.05, edgecolor='k', bottom=df1["contacts"].values)
		elif f == "m1sn.dat":
			df3 = pd.read_csv(f, sep='\t', engine='python', names=["ms1", "ms2", "mm"], skiprows=1)
			s1m_c = []
			for i in range(len(df1["frac"].values)):
				if df1["frac"].values[i] == 1.0:
					s1m_c.append (0)
				else:
					s1m_c.append(df3["ms1"][i]*32/(34*34*34-32-np.floor(df1["frac"][i]*(34*34*34-32))))
			ax.bar(df1["frac"].values, s1m_c, color='royalblue', width=0.05, edgecolor='k', bottom=df1["contacts"].values+df2["s1s2"].values)



	ax.legend(["$s_1$-$s_1$ contacts", "$s_1$-$s_2$ contacts", "$s_1$-$m_1$ contacts"], ncol=3)
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
	ax.axhline (y=0, c='k', linewidth=0.2)
	# ax.yaxis.get_major_locator().set_params(integer=True)
	ax.minorticks_on()
	ax.set_xlim((-0.05, 1.05))
	ax.set_ylim((-1 , 30))
	ax.set_xticks (np.linspace (0, 1, 6))
	ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0f}'))
	plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
	plt.savefig   (args.pn + ".png", bbox_inches='tight', dpi=1200)
