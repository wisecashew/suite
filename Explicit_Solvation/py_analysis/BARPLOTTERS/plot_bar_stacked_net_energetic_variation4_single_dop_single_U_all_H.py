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
parser.add_argument('-H', dest='H', action='store', nargs='+', type=str, help='Provide energy surface.') 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()


if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = aux.dir2U (os.listdir("."))
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.6,2.0), constrained_layout=True)
	rc('font', weight='bold')

	PLOT_DICT  = {}
	ERROR_DICT = {}
	dop            = args.dop
	dump_file      = args.e
	# starting_index = args.s

	ms_max = 25*2+(args.dop-2)*24
	mma_mean = np.asarray([])
	mmn_mean = np.asarray([])
	ms1a_mean = np.asarray([])
	ms1n_mean = np.asarray([])
	ms2a_mean = np.asarray([])
	ms2n_mean = np.asarray([])
	frac_list = []
	H = args.H[0]
	for U in U_list:
		frac_list.append (1- aux.get_frac (U + "/geom_and_esurf.txt" ) )
		mma_list = np.asarray([])
		mmn_list = np.asarray([])
		ms1a_list = np.asarray([])
		ms1n_list = np.asarray([])
		ms2a_list = np.asarray([])
		ms2n_list = np.asarray([])

		num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/E_"+str(args.H[0]) ) ) )

		for num in num_list:
			# print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
			df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/E_"+str(args.H[0])+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
			mma_list  = np.hstack  ( (mma_list, np.mean(df["mm_aligned"  ].values[-2000:]*2/args.dop) ) )
			mmn_list  = np.hstack  ( (mmn_list, np.mean(df["mm_naligned" ].values[-2000:]*2/args.dop) ) )
			ms1a_list  = np.hstack ( (ms1a_list, np.mean(df["ms1_aligned" ].values[-2000:]/args.dop  ) ) )
			ms1n_list  = np.hstack ( (ms1n_list, np.mean(df["ms1_naligned"].values[-2000:]/args.dop  ) ) )
			ms2a_list  = np.hstack ( (ms2a_list, np.mean(df["ms2_aligned" ].values[-2000:]/args.dop  ) ) )
			ms2n_list  = np.hstack ( (ms2n_list, np.mean(df["ms2_naligned"].values[-2000:]/args.dop  ) ) )

		ms1a_mean = np.hstack ( (ms1a_mean, np.mean(ms1a_list) ) )
		ms1n_mean = np.hstack ( (ms1n_mean, np.mean(ms1n_list) ) )
		ms2a_mean = np.hstack ( (ms2a_mean, np.mean(ms2a_list) ) )
		ms2n_mean = np.hstack ( (ms2n_mean, np.mean(ms2n_list) ) )
		mma_mean  = np.hstack ( (mma_mean , np.mean(mma_list)  ) )
		mmn_mean  = np.hstack ( (mmn_mean , np.mean(mmn_list)  ) )

	PLOT_DICT [H] = (np.flip(mma_mean/26 + mmn_mean/26), np.flip(ms2a_mean/26 + ms2n_mean/26), np.flip(ms1a_mean/26 + ms1n_mean/26))
	print("done!", flush=True)
	j = 0


	ax.bar ( np.arange(len(frac_list)), PLOT_DICT[H][0], color ='darkred', width=0.8, edgecolor='k')
	ax.bar ( np.arange(len(frac_list)), PLOT_DICT[H][1], bottom = PLOT_DICT[H][0], color='steelblue', width=0.8, edgecolor='k')
	ax.bar ( np.arange(len(frac_list)), PLOT_DICT[H][2], bottom = PLOT_DICT[H][0]+PLOT_DICT[H][1], color='darkgreen', width=0.8, edgecolor='k')


	for j in range(len(args.H)):
		ax.tick_params ( direction='in', bottom=True, top=True, left=True, right=True, which='both')
		ax.tick_params ( axis='x', labelsize=6.0, direction="in", left="off", labelleft="on", pad=3, labelrotation=0 )
		ax.tick_params ( axis='y', labelsize=8, direction="in", left="off", labelleft="on" )
		ax.axhline (y=0, c='k', linewidth=1)
		# ax[j].minorticks_on()
		ax.set_ylim((-0.01 , 1.01))
		ax.set_xlim((-0.5 , len(frac_list)-0.5))
		ax.set_xticks (np.arange(len(frac_list)))
		ax.set_xticklabels ([ "{:0.2f}".format(1-i) for i in frac_list], weight='bold')
		ax.set_yticks (np.arange(0, 1.2, 0.2))
		ax.set_yticklabels (np.arange (0, 1.4, 0.2), weight='bold')
		ax.yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
		ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
		for f in fig.get_axes():
			f.label_outer()
	
	plt.savefig   (args.pn, dpi=1200)

