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
from scipy.stats import norm
from scipy.optimize import curve_fit


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--U', dest='U', action='store', type=str, help='Name of forcefields.', default=[]) 
parser.add_argument('--frac', dest='frac', action='store', type=str, help='Name of temperature.', default=[]) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

	fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
	# get the entire list of potential energy surfaces 
	fig = plt.figure   (constrained_layout=True)
	ax  = plt.axes() 
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
	# ax.tick_params(axis='y', labelsize=30)
	start = time.time()

	U = args.U
	T = args.frac
	
	i=0
	ms_max = 208 # 25*2+(args.dop-2)*24
	print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
	sc_list = np.asarray([])
	
	print(f"@frac={T}")
	num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

	for num in num_list: 
		df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=1)
		sc_list = np.hstack ( (sc_list, df["ms1s2_tot"].values[-args.s:] ) )

	print(f"Mean of SC = {np.mean(sc_list)}", flush=True)
	print(f"Std. dev. of SC = {np.std(sc_list)}", flush=True)
	print("done!", flush=True)

	i=0
	ax.hist(sc_list, bins=50, density=True)
	hist_values, bin_edges = np.histogram(sc_list, bins=50, density=True)

	
	def beta_fit(x, a, b, loc, scale):
		rv = beta(a, b, loc=loc, scale=scale)
		return rv.pdf(x)

	def chi_fit(x, k, loc, scale):
		rv = chi2(k, loc=loc, scale=scale)
		return rv.pdf(x)

	def poisson_fit(x, k, loc):
		rv = poisson(k, loc=loc)
		return rv.pmf(x)

	mu, std = norm.fit(sc_list)
	print(f"mu = {mu}, std = {std}")
	
	alpha, bbeta, loc, scale = beta.fit(sc_list, floc=0, fscale=26*(34**3)-37*(34**2))
	print(f"alpha = {alpha}, beta = {bbeta}, loc = {loc}, scale = {scale}", flush=True)

	x = np.arange(np.min(sc_list), np.max(sc_list))
	rv = beta(alpha, bbeta, loc=loc, scale=scale)
	print(f"mean = {rv.mean()}, std = {rv.std()}")
	_pdf = rv.pdf(x)
	# ax.plot(x, _pdf, c='green', ls='--')

	p = norm.pdf(x, mu, std)
	ax.plot(x, p, c='k', ls='--')

	# ax.set_xlim(0, 26*(34**3)-37*(34**2))
	ax.set_aspect ('auto')
	fig.savefig (args.pn, bbox_inches='tight', dpi=1200)


