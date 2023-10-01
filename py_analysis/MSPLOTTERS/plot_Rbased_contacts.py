#!/home/satyend/.conda/envs/phase/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as colors
import time
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tck
import argparse 
import sys
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux 
import os 
from pathlib import Path

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.', default=32) 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature list.')
parser.add_argument('--xllim', dest='xllim', action='store', type=float, help='Provide a lower x limit.' , default=0.01)
parser.add_argument('--xulim', dest='xulim', action='store', type=float, help='Provide an upper x limit.', default=100 )
parser.add_argument('--yllim', dest='yllim', action='store', type=float, help='Provide a lower y limit.' , default=0)
parser.add_argument('--yulim', dest='yulim', action='store', type=float, help='Provide an upper y limit.', default=1)
parser.add_argument('--key', dest='key', action='store', type=str, help='Provide the key for which contacts to plot.', default="mm_tot")
parser.add_argument('--database', dest='db', action='store', type=str, help='The database where the flory exponents is stored.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--R', dest='R', metavar='RX', action='store', type=str, help='Name of forcefield in database to plot.')
parser.add_argument('--U', dest='U', metavar='UX', action='store', type=str, help='Name of forcefield directory.')
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 
color_dict = {"R0":'#369DE8', "R1": '#1FB967', "R2":'#B9B41F', "R3":'#B91F72'}

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")

	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(1.7,1.7) )
	lsize = 14
	fig.tight_layout()
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	norm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 )
	fdict = {'color': 'black', 'weight': 'normal', 'size': lsize}

	color1 = np.array([131, 159, 192]) / 255.0  # #839FC0 in RGB
	color2 = np.array([137, 245, 162]) / 255.0   # #ED8151 in RGB
	cmap = colors.LinearSegmentedColormap.from_list('custom', [color2, "white", color1])
	cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)
	y = np.array([0,1])

	df = pd.read_csv (args.db, sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
	df = df.loc[df["U"] == args.R]
	temperatures = df["T"].values
	print (df)

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

        # plotted the background


	start = time.time()

	U_list = [args.U]
	
	PLOT_DICT = dict()
	
	i=0
	if args.key == "ms1_tot" or args.key == "ms1_aligned" or args.key == "ms1_naligned":
		ms_max = 25*2+(args.dop-2)*24
	elif args.key == "mm_tot" or args.key == "mm_aligned" or args.key == "mm_naligned":
		ms_max = 208

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		ms_list = np.asarray([])
		ms_err  = np.asarray([])
		ms_mean = np.asarray([])
		temperatures = args.T
		for temp in temperatures:
			skip = 0
			ms_list = np.asarray ([])
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

			for num in num_list: 
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
				ms_list = np.hstack ( (ms_list, np.mean(df[args.key].values[-2000:] ) ) )

			ms_err  = np.hstack ( (ms_err ,  (np.std (ms_list) / np.sqrt(len(num_list) ) ) ) )
			ms_mean = np.hstack ( (ms_mean,  np.mean (ms_list) ) )

		PLOT_DICT[U] = ( ms_mean, ms_err )
		i += 1
		print("done!", flush=True)


	for U in U_list:
		rgba_color = color_dict[args.R]
		plt.errorbar ( np.array(temperatures), PLOT_DICT[U][0] / ms_max, yerr=PLOT_DICT[U][1] / ms_max, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
		plt.plot     ( np.array(temperatures), PLOT_DICT[U][0] / ms_max, linestyle='--', marker='o',  markeredgecolor='k', linewidth=1, color=rgba_color, label="_nolabel_", markersize=8/1.2, zorder=10, clip_on=False)
		i += 1

	yticks = np.arange(0.0, 1.2, 0.2)
	ax.set_yticks ( yticks )
	ax.set_yticklabels (ax.get_yticks(), fontdict=fdict, font=fpath) 
	ax.set_ylim   (args.yllim, args.yulim)
	ax.set_xlim   (args.xllim, args.xulim)
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticklabels ([]) # [0.01, 0.1, 1.0, 10.0, 100.0], fontdict=fdict, font=fpath)
	ax.set_xscale ("log")
	ax.set_xticks (np.logspace(-2, 2, 5))
	ax.set_xticks ( np.hstack((np.arange(0.01,0.1,0.01), np.arange(0.1, 1, 0.1), np.arange(1,10,1), np.arange(10,100,10))), minor=True)
	# ax.yaxis.set_major_formatter(tck.StrMethodFormatter('{x:1.1f}') )
	ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
	ax.set_yticklabels ([])
	ax.set_xticklabels ([])
	ax.set_aspect ('auto')
	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)

