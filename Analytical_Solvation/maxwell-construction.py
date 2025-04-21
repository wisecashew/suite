import mpmath as mp
import pandas as pd
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",      dest='o',     type=str,   action='store', help="Name of image.")
parser.add_argument ("--EMMA",  dest='emma',  type=float, action='store', help="EMMA value.")
parser.add_argument ("--EMMN",  dest='emmn',  type=float, action='store', help="EMMN value.")
parser.add_argument ("--EMSA",  dest='emsa',  type=float, action='store', help="EMSA value.")
parser.add_argument ("--EMSN",  dest='emsn',  type=float, action='store', help="EMSN value.")
parser.add_argument ("--T",     dest='T',     type=float, action='store', help="temperature value.")
parser.add_argument ("--nbins", dest='nbins', type=int,   action='store', help="number of bins in the histogram.", default=100)
args = parser.parse_args()

class Maxwell:
	def __init__(self, thermo):
		self.Nm   = thermo["Nm"]
		self.z    = 26
		self.pw   = 0.25
		self.kb   = 1
		self.emm_a = args.emma 
		self.ems_a = args.emsa
		self.emm_n = args.emmn
		self.ems_n = args.emsn 
	
	def p_chain(self, Nmm):
		# parameters for beta distribution
		a     =  6.004256058771041
		b     = 27.319457320842343
		loc   = 31
		scale = 208 - 31
		rv = beta(a, b, loc=loc, scale=scale)
		return rv.pdf(Nmm)
		
	comb_mm   = lambda self, Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(self.pw, Nmm_a) * mp.power(1-self.pw, Nmm-Nmm_a)
	comb_ms   = lambda self, Nmm, Nms_a: mp.factorial(Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(Nms_tot(Nmm)-Nms_a)) * mp.power(self.pw, Nms_a) * mp.power(1-self.pw, Nms_tot(Nmm)-Nms_a)
	energy    = lambda self, Nmm, Nmm_a, Nms_a, T: mp.exp(-1/(kb * T) * (Nmm_a * self.emm_a + (Nmm - Nmm_a) * self.emm_n + Nms_a * self.ems_a + (Nms_tot(Nmm) - Nms_a) * self.ems_n))
	boltzmann = lambda self, Nmm, Nmm_a, Nms_a, T: p_chain(Nmm) * comb_mm(Nmm, Nmm_a) * comb_ms(Nmm, Nms_a) * mp.exp(-1/T * self.energy(Nmm, Nmm_a, Nms_a))
	


if __name__=="__main__":

	# get some time counters 
	start_global = time.time()

	fig = plt.figure(figsize=(2,2))
	ax  = plt.axes()

	# We change the fontsize of minor ticks label 
	ax.tick_params(axis='both', which='major', labelsize=6)
	ax.tick_params(axis='both', which='major', labelsize=6)
	start = time.time()

	info    = dict()
	Nm      = 32
	z       = 26
	kb      = 1
	pw      = 0.25
	emm_a   = args.emma
	ems_a   = args.emsa
	emm_n   = args.emmn
	ems_n   = args.emsn
	T       = args.T

	energy    = lambda Nmm, Nmm_a, Nms_a: mp.exp(-1/(kb * T) * (Nmm_a * emm_a + (Nmm - Nmm_a) * emm_n + Nms_a * ems_a + (Nms_tot(Nmm) - Nms_a) * ems_n))

	info["T"]  = T
	Nmm_list   = []
	Nms_list   = []

	# parameters for beta distribution
	a     =  6.004256058771041
	b     = 27.319457320842343
	loc   = 31
	scale = 208 - 31

	rv = beta(a, b, loc=loc, scale=scale)

	Nms_tot   = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain   = lambda Nmm: rv.pdf(Nmm)
	comb_mm   = lambda Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(pw, Nmm_a) * mp.power(1-pw, Nmm-Nmm_a)
	comb_ms   = lambda Nmm, Nms_a: mp.factorial(Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(Nms_tot(Nmm)-Nms_a)) * mp.power(pw, Nms_a) * mp.power(1-pw, Nms_tot(Nmm)-Nms_a)
	energy    = lambda Nmm, Nmm_a, Nms_a: (Nmm_a * emm_a + (Nmm - Nmm_a) * emm_n + Nms_a * ems_a + (Nms_tot(Nmm) - Nms_a) * ems_n)
	boltzmann = lambda Nmm, Nmm_a, Nms_a: p_chain(Nmm) * comb_mm(Nmm, Nmm_a) * comb_ms(Nmm, Nms_a) * mp.exp(-1/T * energy(Nmm, Nmm_a, Nms_a))

	energy_store  = []
	weights_store = []
	for Nmm in range(31, 209, 10):
		start = time.time()
		print(f"\t@ Nmm = {Nmm}...", flush=True, end=' ')

		Nmm_a_vals = range(0, Nmm)
		Nms_a_vals = range(0, Nms_tot(Nmm))

		# use list comprehension to avoid overheads

		energies = [energy(Nmm, Nmm_a, Nms_a) for Nmm_a in Nmm_a_vals for Nms_a in Nms_a_vals]
		weights  = [boltzmann(Nmm, Nmm_a, Nms_a) for Nmm_a in Nmm_a_vals for Nms_a in Nms_a_vals]

		energy_store.extend(energies)
		weights_store.extend(weights)
		stop = time.time()
		print(f"time: {stop - start} seconds.", flush=True)

	print(f"The minimum in energy is {np.min(energy_store)}", flush=True)
	print(f"The maximum in energy is {np.max(energy_store)}", flush=True)

	Emin = np.min(energy_store)
	Emax = np.max(energy_store)

	n_bins = args.nbins
	bins   = np.linspace(Emin, Emax, n_bins+1)
	hist, bin_edges = np.histogram(energy_store, bins=bins, weights=weights_store, density=True)

	ax.set_xlim(Emin, Emax)
	ax.minorticks_on()
	ax.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor='k', align='edge')
	ax.set_xlabel('Energy', fontsize=9)
	ax.set_ylabel('Total Weight', fontsize=9)
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
	stop_global = time.time()

	print(f"Time for computation is {stop_global-start_global} seconds.", flush=True)
