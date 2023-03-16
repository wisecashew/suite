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
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.')

args = parser.parse_args()


if __name__=="__main__":

	# get the entire list of potential energy surfaces 
	U_list = args.U
	fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4.6,2.0), constrained_layout=True)
	rc ('font', weight='bold')

	PLOT_DICT  = {}
	ERROR_DICT = {}
	dop            = args.dop
	dump_file      = args.e
	# starting_index = args.s

	ms_max = 25*2+(args.dop-2)*24
	for U in U_list:

		temperatures = [0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0]
		temperatures.sort()
		mma_list = np.asarray([])
		mma_mean = np.asarray([])
		mmn_list = np.asarray([])
		mmn_mean = np.asarray([])
		msa_list = np.asarray([])
		msa_mean = np.asarray([])
		msn_list = np.asarray([])
		msn_mean = np.asarray([])

		for T in temperatures: 

			mma_list = np.asarray ([]) 
			mmn_list = np.asarray ([])
			msa_list  = np.asarray ([]) 
			msn_list  = np.asarray ([]) 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

			for num in num_list:
				# print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
				mma_list  = np.hstack ( (mma_list, np.mean(df["mm_aligned"  ].values[-2000:]*2/args.dop) ) )
				mmn_list  = np.hstack ( (mmn_list, np.mean(df["mm_naligned" ].values[-2000:]*2/args.dop) ) )
				msa_list  = np.hstack ( (msa_list, np.mean(df["ms1_aligned" ].values[-2000:]/args.dop  ) ) )
				msn_list  = np.hstack ( (msn_list, np.mean(df["ms1_naligned"].values[-2000:]/args.dop  ) ) )

			msa_mean = np.hstack ( (msa_mean, np.mean(msa_list) ) )
			msn_mean = np.hstack ( (msn_mean, np.mean(msn_list) ) )
			mma_mean = np.hstack ( (mma_mean, np.mean(mma_list) ) )
			mmn_mean = np.hstack ( (mmn_mean, np.mean(mmn_list) ) )

		PLOT_DICT [U] = (mma_mean/26, mmn_mean/26, msa_mean/26, msn_mean/26)
		print("done!", flush=True)
	j = 0
	for U in U_list:
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][0], color ='darkred', width=0.8, edgecolor='k')
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][1], bottom = PLOT_DICT[U][0], color='lightcoral', width=0.8, edgecolor='k')
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][2], bottom = PLOT_DICT[U][0]+PLOT_DICT[U][1], color='steelblue', width=0.8, edgecolor='k')
		ax[j].bar ( np.arange(len(temperatures)), PLOT_DICT[U][3], bottom = PLOT_DICT[U][0]+PLOT_DICT[U][1]+PLOT_DICT[U][2], color='lightskyblue', width=0.8, edgecolor='k')
		j += 1
	###
	# write data points out for text later 
	g = open ("m1sn.dat", 'w')
	g.write ("m-m-a	m-m-n	m-s-a	m-s-n\n")
	for i in range(len (PLOT_DICT[U][0])):
		g.write ("{}	{}	{}	{}\n".format(PLOT_DICT[U][0][i], PLOT_DICT[U][1][i], PLOT_DICT[U][2][i], PLOT_DICT[U][3][i] ) )
	g.close()
	for j in range(len(U_list)):
	# 	if j == 0:
	# 		ax[j].text (x=-0.5, y=1.03, s="$\eta ^{a}$=-0.1\n$\eta^{a'}$=-0.1", fontsize=7)
	# 		ax[j].set_ylabel ("$P(i|m)$", fontsize=12)
	# 	else:
	# 		ax[j].text (x=-0.5, y=1.03, s="$\eta ^{a}$=0.1\n$\eta^{a'}$=-0.1", fontsize=7)
	# 	if j == 1:
	# 		ax[j].legend(["$i \\rightarrow m$ (aligned)", "$i \\rightarrow m$ (misaligned)", "$i \\rightarrow s$ (aligned)", "$i \\rightarrow s$ (misaligned)"], ncol=1, loc="upper left", bbox_to_anchor=(1,1), frameon=False, fontsize=4) 
		ax[j].tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
		ax[j].tick_params ( axis='x', labelsize=6.0, direction="in", left="off", labelleft="on", pad=3, labelrotation=45 )
		ax[j].tick_params ( axis='y', labelsize=8.0, direction="in", left="off", labelleft="on" )
		ax[j].axhline (y=0, c='k', linewidth=1)
		# ax[j].minorticks_on()
		ax[j].set_ylim((-0.01 , 1.01))
		ax[j].set_xlim((-0.5, len(temperatures)-0.5))
		ax[j].set_xticks (np.arange(len(temperatures)))
		ax[j].set_xticklabels ([str(i) for i in temperatures], weight='bold')
		ax[j].set_yticks (np.arange(0, 1.2, 0.2))
		ax[j].set_yticklabels (np.arange (0, 1.4, 0.2), weight='bold')
		ax[j].yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
		ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator())
		# ax[j].set_xlabel ("$T^*$", fontsize=12)
		for f in fig.get_axes():
			f.label_outer()
	
	plt.savefig   (args.pn, dpi=1200)

