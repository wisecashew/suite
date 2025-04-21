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
parser.add_argument('--dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('--dop-suffix', dest='dsuf', action='store', type=str, help='Provide suffix of dop file.', default="") 
parser.add_argument('--U', dest='U', action='store', nargs='+', type=str, help='Name of forcefields.', default=[]) 
parser.add_argument('--solv', dest='solv', metavar='solvation', action='store', type=str, help='Name of energy file to parse information.', default='energy') 
parser.add_argument('--energy', dest='energy', metavar='energy', action='store', type=str, help='Name of solvation file to parse information.', default='solvation') 
parser.add_argument('--xllim',   dest='xllim', action='store', type=float, help="Lower x limit.", default=None)
parser.add_argument('--xulim',   dest='xulim', action='store', type=float, help="Upper x limit.", default=None)
parser.add_argument('--yllim',   dest='yllim', action='store', type=float, help="Lower y limit.", default=None)
parser.add_argument('--yulim',   dest='yulim', action='store', type=float, help="Upper y limit.", default=None)
parser.add_argument('--color', dest='col', metavar='color', nargs='+', action='store', type=str, help='Name of color of file.', default='coral') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
	# get the entire list of potential energy surfaces 
	fig = plt.figure(figsize=(1.65,1.65))
	ax  = plt.axes()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='y', labelsize=6)
	ax.tick_params(axis='x', labelsize=6)
	start = time.time()

	U_list = args.U
	
	i = 0
	Z = 26
	scale = args.dop * Z
	
	solv_list = np.asarray([])
	ms_list   = np.asarray([])

	for U in U_list:
		print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
		x_vals    = aux.dir2float(os.listdir(str(U)+"/DOP_"+str(args.dop)+args.dsuf))
		for x in x_vals:
			skip = 0
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop) + args.dsuf + "/" + str(x) ) ) )

			for num in num_list:
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop) + args.dsuf + "/"+str(x)+"/"+args.solv+"_"+str(num)+".mc", sep=' \| ', names=["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "timestep"], engine='python')
				solv_list = np.hstack ((solv_list, df["total"].values))
				df = pd.read_csv(str(U)+"/DOP_"+str(args.dop) + args.dsuf + "/"+str(x)+"/"+args.energy+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python')
				ms_list = np.hstack((ms_list, df["ms1_tot"].values))

		i += 1
		print("done!", flush=True)

	ax.scatter(ms_list, solv_list, marker='o',  ec='k', s=8/1.3, zorder=10, clip_on=False, c=args.col[0])

	ax.minorticks_on()
	fig.savefig (args.pn, bbox_inches='tight', dpi=1200)


