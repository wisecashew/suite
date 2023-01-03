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
from matplotlib import rc,rcParams

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-U', dest='U', action='store', nargs='+', type=str, help='Provide energy surface.') 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()


if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = args.U
	fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(1.8,1.5)) # , constrained_layout=True)
	fig.subplots_adjust (hspace=0.3)
	rc('font', weight='bold')
	# ax[0].get_shared_x_axes().join(ax[0], ax[1])
	# instantiate plt figure 
	# fig.figure( figsize=(8,6) )

	PLOT_DICT  = {}
	ERROR_DICT = {}
	dop            = args.dop
	dump_file      = args.e
	# starting_index = args.s

	ms_max = 25*2+(args.dop-2)*24
	for U in U_list:
        
        # temperatures = aux.dir2float (os.listdir(str(U)+"/DOP_32"))
		temperatures = [0.01, 0.1, 0.5, 1.0, 2.5, 10.0, 25.0, 50.0, 100.0]
		temperatures.sort()
		mm_list = np.asarray([])
		mm_mean = np.asarray([])
		ms_list = np.asarray([])
		ms_mean = np.asarray([])

		for T in temperatures: 
			mm_list = np.asarray ([]) 
			ms_list  = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

			for num in num_list: 
				# print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
				mm_list   = np.hstack ( (mm_list, np.mean(df["mm_aligned"  ].values[-2000:]*2/args.dop + df["mm_naligned"].values[-2000:]*2/args.dop ) ) )
				ms_list   = np.hstack ( (ms_list, np.mean(df["ms1_aligned" ].values[-2000:]/args.dop + df["ms1_naligned" ].values[-2000:]/args.dop ) ) )

			ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )
			mm_mean = np.hstack ( (mm_mean, np.mean(mm_list) ) )

		PLOT_DICT [U] = (mm_mean, ms_mean)
		print("done!", flush=True)
    
	j = 0
	for U in U_list:
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][0]/26, color ='darkred', width=0.8, edgecolor='k')
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][1]/26, bottom = PLOT_DICT[U][0]/26, color='steelblue', width=0.8, edgecolor='k')
		j += 1
	###
	# write data points out for text later 
	g = open ("m1sn.dat", 'w')
	g.write ("m-m		m-s\n")
	for i in range(len (PLOT_DICT[U][0])):
		g.write ("{}	{}\n".format(PLOT_DICT[U][0][i], PLOT_DICT[U][1][i] ) )
	g.close()
	for j in range(len(U_list)):
		if j == 0:
			ax[j].legend(["$i \\rightarrow m$", "$i \\rightarrow s$"], ncol=2, bbox_to_anchor=(0.45, 1.25), loc="upper left", frameon=False, fontsize=4) # , bbox_to_anchor=(0.5, 1.13))
		ax[j].tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both' )
		ax[j].tick_params ( axis='x', labelsize=4, direction="in", left="off", labelleft="on", labelrotation=45, pad=1 )
		ax[j].tick_params ( axis='y', labelsize=4, direction="in", left="off", labelleft="on" )
		ax[j].axhline (y=0, c='k', linewidth=1)
		ax[j].minorticks_on()
		if j == 0:
			ax[j].set_ylim((0 , 1))
		else:
			ax[j].set_ylim((0 , 1))
		ax[j].set_xlim((-0.5, len(temperatures)-0.5))
		ax[j].set_xticks (np.arange(len(temperatures)))
		ax[j].set_xticklabels ([str(i) for i in temperatures], weight='bold')
		ax[j].set_yticks (np.arange (0, 1.2, 0.2))
		ax[j].set_yticklabels (np.arange (0, 1.2, 0.2), weight='bold')
		ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator())
		ax[j].yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
		for f in fig.get_axes():
			f.label_outer()
	plt.savefig   (args.pn, dpi=1200)

