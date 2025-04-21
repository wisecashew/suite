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
from scipy.stats import beta
from scipy.stats import chi2
from scipy.stats import poisson
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
	# get the entire list of potential energy surfaces 
	fig = plt.figure   ( figsize=(1.7,1.7), constrained_layout=True )
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax.tick_params(axis='y', labelsize=4)
	ax.tick_params(axis='x', labelsize=4)
	start = time.time()

	i=0
	dop = args.dop
	mm_list = np.asarray([])
	
	df = pd.read_csv(args.e, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", \
	"ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=1)
	mm_list = np.hstack ( (mm_list, df["mm_tot"].values ) )

	print(f"Mean of MM = {np.mean(mm_list)}", flush=True)
	print("done!", flush=True)

	i=0
	ax.hist(mm_list, bins=50, density=True)
	hist_values, bin_edges = np.histogram(mm_list, bins=50, density=True)

	def beta_fit(x, a, b, loc, scale):
		rv = beta(a, b, loc=loc, scale=scale)
		return rv.pdf(x)

	def chi_fit(x, k, loc, scale):
		rv = chi2(k, loc=loc, scale=scale)
		return rv.pdf(x)

	def poisson_fit(x, k, loc):
		rv = poisson(k, loc=loc)
		return rv.pmf(x)

	print(np.max(mm_list), np.min(mm_list))
	alpha, bbeta, loc, scale = beta.fit(mm_list, floc=0, fscale=5*dop)
	print(f"beta fit: alpha = {alpha} , beta = {bbeta},   loc = {loc},   scale = {scale}",   flush=True)

	x = np.arange(np.min(mm_list), np.max(mm_list))
	rv = beta(alpha, bbeta, loc=loc, scale=scale)
	_pdf = rv.pdf(x)
	ax.plot(x, _pdf, c='green', ls='-.', label="fixed loc, scale")

	ax.set_xlim(0, np.max(mm_list))
	ax.legend(loc="upper right", fontsize=4)
	ax.set_aspect ('auto')
	fig.savefig (args.pn, bbox_inches='tight', dpi=1200)

	scaled_mm_list = mm_list/((2.19*dop-7.49)*10/3)
	salpha, sbeta, sloc, sscale = beta.fit(scaled_mm_list, floc=0, fscale=1)
	print(f"beta fit: alpha = {salpha} , beta = {sbeta},   loc = {sloc},   scale = {sscale}",   flush=True)

	fig0 = plt.figure   (num=2, figsize=(1.7,1.7), constrained_layout=True )
	ax0  = plt.axes()
	ax0.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	ax0.tick_params(axis='y', labelsize=4)
	ax0.tick_params(axis='x', labelsize=4)

	x = np.linspace(0, 1, 1000)
	rv = beta(14.819, 34.073599, loc=sloc, scale=sscale)
	_pdf = rv.pdf(x)
	print(f"mean for between 0 and 1 = {rv.mean()}")
	ax0.hist(scaled_mm_list, bins=50, density=True)
	ax0.plot(x, _pdf, c='green', ls='-.', label="fixed loc, scale")
	ax0.set_xlim(0, 1)
	ax0.legend(loc="upper right", fontsize=4)
	ax0.set_aspect ('auto')
	fig0.savefig (args.pn+"_scaled", bbox_inches='tight', dpi=1200)
