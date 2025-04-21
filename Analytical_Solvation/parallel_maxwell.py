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
parser.add_argument ("--nproc", dest='nproc', type=int, action='store', help="number of parallel processors.")
parser.add_argument ("--nbins", dest='nbins', type=int,   action='store', help="number of bins in the histogram.", default=100)
args = parser.parse_args()

class Maxwell:
	def __init__(self, Nm, T, nproc):
		self.Nm   = Nm
		self.z    = 26
		self.pw   = 0.25
		self.kb   = 1
		self.emm_a = args.emma 
		self.ems_a = args.emsa
		self.emm_n = args.emmn
		self.ems_n = args.emsn 
		self.T     = T
		self.nproc = nproc
	
	def p_chain(self, Nmm):
		# parameters for beta distribution
		a     =  6.004256058771041
		b     = 27.319457320842343
		loc   = 31
		scale = 208 - 31
		rv = beta(a, b, loc=loc, scale=scale)
		return rv.pdf(Nmm)
		
	Nms_tot   = lambda self, Nmm: 770 - 2 * (Nmm - (self.Nm - 1))
	comb_mm   = lambda self, Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(self.pw, Nmm_a) * mp.power(1-self.pw, Nmm-Nmm_a)
	comb_ms   = lambda self, Nmm, Nms_a: mp.factorial(self.Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(self.Nms_tot(Nmm)-Nms_a)) * mp.power(self.pw, Nms_a) * mp.power(1-self.pw, self.Nms_tot(Nmm)-Nms_a)
	energy    = lambda self, Nmm, Nmm_a, Nms_a: (Nmm_a * self.emm_a + (Nmm - Nmm_a) * self.emm_n + Nms_a * self.ems_a + (self.Nms_tot(Nmm) - Nms_a) * self.ems_n)
	boltzmann = lambda self, Nmm, Nmm_a, Nms_a: self.p_chain(Nmm) * self.comb_mm(Nmm, Nmm_a) * self.comb_ms(Nmm, Nms_a) * mp.exp(-1/self.T * self.energy(Nmm, Nmm_a, Nms_a))
	
	def run_partial_sum(self, Nmm_range):
		energy_store  = []
		weights_store = []

		for Nmm in Nmm_range:
			print(f"@ Nmm = {Nmm}...", flush=True)
			Nmm_a_vals = range(0, Nmm)
			Nms_a_vals = range(0, self.Nms_tot(Nmm))

			# use list comprehension to avoid overheads

			energies = [self.energy(Nmm, Nmm_a, Nms_a) for Nmm_a in Nmm_a_vals for Nms_a in Nms_a_vals]
			weights  = [self.boltzmann(Nmm, Nmm_a, Nms_a) for Nmm_a in Nmm_a_vals for Nms_a in Nms_a_vals]

			energy_store.extend(energies)
			weights_store.extend(weights)
		
		return (energy_store, weights_store)

	def run_total_sum(self):

		# get the ranges
		total_range = np.arange(31, 209)
		seg_length  = len(total_range) // self.nproc

		# set up some preprocessing
		if len(total_range) % self.nproc == 0:
			pass
		else:
			seg_length += 1
		
		# get up more containers
		segments    = list()
		for i in range(self.nproc):
			segments.append(total_range[i::self.nproc])
		
		print(f"Length of segments is {len(segments)}", flush=True)
		# make the pool
		pool = multiprocessing.Pool(processes=self.nproc)
		
		# run the parallel calc
		results = pool.starmap(self.run_partial_sum, zip(segments)) 

		# set up the containers
		total_energy_store = list()
		total_weight_store = list()

		# grab all the results
		for res in results:
			total_energy_store.extend(res[0])
			total_weight_store.extend(res[1])

		return total_energy_store, total_weight_store

if __name__=="__main__":

	# get some time counters 
	start_global = time.time()

	# define properties of objects
	Nm      = 32
	T       = args.T
	nproc   = args.nproc
	maxwell = Maxwell(Nm, T, nproc)

	# run the big summation	
	energy_store, weights_store = maxwell.run_total_sum()

	# print out the relevant stuff
	print(f"The minimum in energy is {np.min(energy_store)}", flush=True)
	print(f"The maximum in energy is {np.max(energy_store)}", flush=True)

	# get the lowest and highest energies
	Emin = np.min(energy_store)
	Emax = np.max(energy_store)

	n_bins = args.nbins
	bins   = np.linspace(Emin, Emax, n_bins+1)
	hist, bin_edges = np.histogram(energy_store, bins=bins, weights=weights_store, density=True)

	# set up the plot 
	fig = plt.figure(figsize=(2, 2))
	ax  = plt.axes()

	# We change the fontsize of minor ticks label 
	ax.tick_params(axis='both', which='major', labelsize=6)
	ax.tick_params(axis='both', which='major', labelsize=6)
	ax.set_xlim(Emin, Emax)
	ax.minorticks_on()
	ax.bar(bin_edges[:-1], hist, edgecolor='k', align='edge', width=1)
	ax.set_xlabel('Energy', fontsize=9)
	ax.set_ylabel('Total Weight', fontsize=9)
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")

	# get the timings
	stop_global = time.time()
	print(f"Time for computation is {stop_global-start_global} seconds.", flush=True)
