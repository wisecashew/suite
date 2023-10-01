#!/home/satyend/.conda/envs/phase/bin/python

import pandas as pd 
from pathlib import Path
import numpy as np 
import matplotlib
from matplotlib.ticker import StrMethodFormatter
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as tck
import argparse 
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux 
import os 
from matplotlib import rc,rcParams


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int,  help='Provide degree of polymerization.')
parser.add_argument('--database', dest='db', action='store', type=str,  help='Provide location of database.')
parser.add_argument('--figsize', dest='figsize', action='store', nargs=2, type=float, help='Enter dimensions of image (width, height) (default: (2.5, 2.5)', default=[2.5,2.5])
parser.add_argument('--yulim', dest='yulim', action='store', type=float, help='Provide an upper limit to y-axis.', default=None)
parser.add_argument('--yllim', dest='yllim', action='store', type=float, help='Provide a lower limit to y-axis.',  default=None)
parser.add_argument('--ylabels', dest='ylabels', action='store_true', default=False, help='Provide this label if you want to see ylabels.')
parser.add_argument('--U',  dest='U',   action='store', nargs='+', type=str,       help='Provide energy surface.')
parser.add_argument('--T',  dest='T',   action='store', nargs='+', type=float,     help='Provide temperatures to probe.')
parser.add_argument('--R', dest='R', metavar='RX', action='store', type=str, help='Name of forcefield in database to plot.')
parser.add_argument('--key', dest='key', action='store', nargs='+', type=str, help='Enter what kind of plot you want. \
\nOptions: 1. mm_contacts_over_total, \n2. ms_contacts_over_total, \n3. mm_aligned_fraction, \n4. ms_aligned_fraction, \n5. aligned_over_total, \n6. mm_over_ms, 7. naligned_over_aligned.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()


if __name__=="__main__":

	fsize = args.figsize
	color_dict = {"R0":'#369DE8', "R1": '#1FB967', "R2":'#B9B41F', "R3":'#B91F72'}

	key = args.key

	if key is None:
		print ("Please provide a key.")
		exit ()

	fig = plt.figure(figsize=(fsize[0],fsize[1]))
	ax  = plt.axes ()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5)



	U_list = args.U

	PLOT_DICT  = {}
	ERROR_DICT = {}
	dop            = args.dop
	dump_file      = args.e
	# starting_index = args.s

	ms_max = 25*2+(args.dop-2)*24
	for U in U_list:

		temperatures = args.T
		mma_list = np.array([])
		mma_mean = np.array([])
		mmn_list = np.array([])
		mmn_mean = np.array([])
		msa_list = np.array([])
		msa_mean = np.array([])
		msn_list = np.array([])
		msn_mean = np.array([])

		for T in temperatures: 

			mma_list = np.array ([])
			mmn_list = np.array ([])
			msa_list = np.array ([])
			msn_list = np.array ([])
			num_list = np.unique (aux.dir2nsim(os.listdir(str(U)+"/DOP_"+str(args.dop)+"/"+str(T))))

			for num in num_list:
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
				mma_list  = np.hstack ((mma_list, np.mean(df["mm_aligned"  ].values[-2000:])))
				mmn_list  = np.hstack ((mmn_list, np.mean(df["mm_naligned" ].values[-2000:])))
				msa_list  = np.hstack ((msa_list, np.mean(df["ms1_aligned" ].values[-2000:])))
				msn_list  = np.hstack ((msn_list, np.mean(df["ms1_naligned"].values[-2000:])))

			msa_mean = np.hstack ( (msa_mean, np.mean(msa_list) ) )
			msn_mean = np.hstack ( (msn_mean, np.mean(msn_list) ) )
			mma_mean = np.hstack ( (mma_mean, np.mean(mma_list) ) )
			mmn_mean = np.hstack ( (mmn_mean, np.mean(mmn_list) ) )

		PLOT_DICT [U] = (mma_mean, mmn_mean, msa_mean, msn_mean)
		print("done!", flush=True)
	j = 0
	for U in U_list:
		total = PLOT_DICT[U][0] + PLOT_DICT[U][1] + PLOT_DICT[U][2] + PLOT_DICT[U][3]
		mm    = PLOT_DICT[U][0] + PLOT_DICT[U][1]
		mma   = PLOT_DICT[U][0]
		ms    = PLOT_DICT[U][2] + PLOT_DICT[U][3]
		msa   = PLOT_DICT[U][2]
		ali   = PLOT_DICT[U][0] + PLOT_DICT[U][2]
		nali  = PLOT_DICT[U][1] + PLOT_DICT[U][3]
		rgba_color = color_dict[args.R]
		for key in args.key:
			if key == "mm_contacts_over_total":
				Y = mm/total
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "ms_contacts_over_total":
				Y = ms/total
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "mm_aligned_fraction":
				Y = mma/mm
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "ms_aligned_fraction":
				Y = msa/ms
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "aligned_over_total":
				Y = ali/total
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "naligned_over_total":
				Y = nali/total
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "mm_over_ms":
				Y = mm/ms
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=False, c=rgba_color, zorder=15, ls='--')

			elif key == "naligned_over_aligned":
				Y = nali/ali
				ax.plot(temperatures, Y, marker='o', mec='k', markersize=8/1.3, clip_on=True, c=rgba_color, zorder=15, ls='--')

			else:
				print ("Bad key.")
				exit ()

	# plot the background
	color1 = np.array([131, 159, 192]) / 255.0   # #839FC0 in RGB
	color2 = np.array([137, 245, 162]) / 255.0   # #ED8151 in RGB
	cmap = colors.LinearSegmentedColormap.from_list('custom', [color2, "white", color1])
	cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)

	if args.yulim == None:
		# yulim = np.max(Y) * 1.1 if np.max(Y) > 0 else np.max(Y)*0.9
		yulim = np.max(Y) if np.max(Y) > 0 else np.max(Y)*0.9
	else:
		yulim = args.yulim

	if args.yllim == None:
		# yllim = np.min(Y) * 0.9 if np.max(Y) > 0 else np.min(Y) * 1.1
		yllim = np.min(Y) if np.max(Y) > 0 else np.min(Y) * 1.1
	else:
		yllim = args.yllim

	y = np.array([yllim, yulim])
	
	df = pd.read_csv (args.db, sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
	df = df.loc[df["U"] == args.R]
	temperatures = df["T"].values

	x_old = temperatures
	y_old = df["nu_mean"].values/2

	x_pred = np.logspace (-2, 2, 10000)
	y_pred = np.interp (x_pred, x_old, y_old)

	col_dict = dict()

	for i in range (len(x_pred)):
		col_dict[ x_pred[i] ] = y_pred[i]

	X, Y = np.meshgrid (x_pred, y)
	Z = np.zeros ((2, len(x_pred)))

	for i in range (len(x_pred)):
		for j in range (2):
			Z[j, i]= col_dict [x_pred[i]]

	ax.pcolormesh (X, Y, Z, cmap=cmap, norm=cnorm, shading="auto")


	###
	# write data points out for text later 
	# yticks = np.arange(0.0, 1.1, 0.1) 
	# ax.set_yticks ( yticks )
	if args.ylabels:
		pass
	else:
		ax.set_yticklabels([])
	ax.set_ylim   (yllim, yulim)
	ax.set_xlim   ( 0.01, 100 )
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xscale('log')
	ax.set_xticks (np.logspace(-2,2,5))
	ax.set_xticks ( np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	ax.set_xticklabels([])

	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	ax.set_aspect('auto')

	plt.savefig   (args.pn, dpi=1200, bbox_inches="tight")

