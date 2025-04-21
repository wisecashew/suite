#!/home/satyend/.conda/envs/phase/bin/python

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

parser = argparse.ArgumentParser(description="Calculating the single chain partition function.")
parser.add_argument ("-o",      dest='o',     type=str,   action='store', help="Name of image.")
parser.add_argument ("--EMMA",  dest='emma',  type=float, action='store', help="EMMA value.")
parser.add_argument ("--EMMN",  dest='emmn',  type=float, action='store', help="EMMN value.")
parser.add_argument ("--EMSA",  dest='emsa',  type=float, action='store', help="EMSA value.")
parser.add_argument ("--EMSN",  dest='emsn',  type=float, action='store', help="EMSN value.")
parser.add_argument ("--T",     dest='T',     type=float, action='store', help="temperature value.", nargs='+')
parser.add_argument ("--nproc", dest='nproc', type=int,   action='store', help="number of parallel processors.")
parser.add_argument ("--nbins", dest='nbins', type=int,   action='store', help="number of bins in the histogram.", default=100)
parser.add_argument ("--style", dest='style', type=str,   action='store', help="Enter the computation you want to run (available options: 'distribution', 'contacts').", default=None)
parser.add_argument ("--csv",   dest='csv',   type=str,   action='store', help="Enter name of csv file to dump information into.", default=None)
args = parser.parse_args()

class Maxwell:
	def __init__(self, Nm, T_range, nproc):
		self.Nm      = Nm
		self.z       = 26
		self.pw      = 0.25
		self.kb      = 1
		self.emm_a   = args.emma 
		self.ems_a   = args.emsa
		self.emm_n   = args.emmn
		self.ems_n   = args.emsn 
		self.T       = T_range[0]
		self.T_range = T_range
		self.nproc   = nproc
	
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
	
	# this is to get energy distributions
	def run_partial_sum(self, Nmm_range):
		energy_store  = []
		weights_store = []

		for Nmm in Nmm_range:
			print(f"\t@ Nmm = {Nmm}...", flush=True)
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

	def plot_distributions(self):
		idx = 0
		for T in self.T_range:
			print(f"@T = {T}...", flush=True)
			self.T = T
			energy_store, weights_store = self.run_total_sum()

			# get the lowest and highest energies
			Emin = np.min(energy_store)
			Emax = np.max(energy_store)

			n_bins = args.nbins
			bins   = np.linspace(Emin, Emax, n_bins+1)
			hist, bin_edges = np.histogram(energy_store, bins=bins, weights=weights_store, density=True)

			# set up the plot 
			fig = plt.figure(num=idx, figsize=(2, 2))
			ax  = plt.axes()

			# We change the fontsize of minor ticks label 
			ax.tick_params(axis='both', which='major', labelsize=6)
			ax.tick_params(axis='both', which='major', labelsize=6)
			ax.set_xlim(Emin, Emax)
			ax.minorticks_on()
			ax.bar(bin_edges[:-1], hist, edgecolor='k', align='edge', width=1)
			ax.set_xlabel('Energy', fontsize=9)
			ax.set_ylabel('Total Weight', fontsize=9)
			fig.savefig(args.o+f"_{T}.png", dpi=1200, bbox_inches="tight")

		return 

	# this is to get heat capacity
	def plot_variance(self):
		Cv_store = list()

		for T in self.T_range:
			print(f"@T = {T}", flush=True)
			self.T = T 
			energy_store, weights_store = self.run_total_sum()
			total_weight = mp.fsum(weights_store)
			weights_store = [x/total_weight for x in weights_store]
			aveE = np.sum(np.array(energy_store) * np.array(weights_store))/np.sum(np.array(weights_store))
			Cv   = np.sum(np.array(weights_store) * (np.array(energy_store) - aveE) ** 2)/np.sum(np.array(weights_store))
			Cv_store.append(Cv/(self.kb*self.T**2))

		# set up the plot 
		fig = plt.figure(figsize=(2, 2))
		ax  = plt.axes()

		# We change the fontsize of minor ticks label 
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.plot(self.T_range, Cv_store, marker='o', mec='k', c='lavender', markersize=6, lw=2, ls='--')
		ax.minorticks_on()
		fig.savefig(args.o, dpi=1200, bbox_inches="tight")

		f = open(args.csv, 'w')
		f.write(f"Temperature | Heat_capacity\n")
		for i in range(len(self.T_range)):
			f.write(f"{self.T_range[i]} | {Cv_store[i]}\n")
		f.close()

		return

	def plot_energy(self):
		E_store = list()
		for T in self.T_range:
			print(f"@T = {T}...", flush=True)
			self.T = T 
			energy_store, weights_store = self.run_total_sum()
			total_weight = mp.fsum(weights_store)
			weights_store = [x/total_weight for x in weights_store]
			aveE = np.sum(np.array(energy_store) * np.array(weights_store))/np.sum(np.array(weights_store))
			E_store.append(aveE)

		# set up the plot
		fig = plt.figure(figsize=(2,2))
		ax  = plt.axes()

		# We change the fontsize of minor ticks label 
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.plot(self.T_range, E_store, marker='o', mec='k', c='silver', markersize=6, lw=2, ls='--')
		ax.minorticks_on()
		fig.savefig(args.o, dpi=1200, bbox_inches="tight")

		f = open(args.csv, 'w')
		f.write(f"Temperature | Energy\n")
		for i in range(len(self.T_range)):
			f.write(f"{self.T_range[i]} | {E_store[i]}\n")
		f.close()

		return

	def plot_free_energy(self):
		free_energy_store = list()
		for T in self.T_range:
			print(f"@T = {T}...", flush=True)
			self.T = T
			energy_store, weights_store = self.run_total_sum()
			total_weight = mp.fsum(weights_store)
			free_energy  = - T * mp.log(total_weight)
			free_energy_store.append(free_energy)

		# set up the plot
		fig = plt.figure(figsize=(2,2))
		ax  = plt.axes()

		# We change the fontsize of minor ticks label
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.plot(self.T_range, free_energy_store, marker='o', mec='k', c='darkorange', markersize=6, lw=2, ls='--')
		ax.minorticks_on()
		fig.savefig(args.o, dpi=1200, bbox_inches="tight")

		f = open(args.csv, 'w')
		f.write(f"Temperature | Free_Energy\n")
		for i in range(len(self.T_range)):
			f.write(f"{self.T_range[i]} | {free_energy_store[i]}\n")
		f.close()

		return	

	# this is to get Nmm contacts
	def run_Nmm_sum(self, Nmm_range):

		ave_Nmm = 0
		Z       = 0

		for Nmm in Nmm_range:
			for Nmm_a in range(0, Nmm):
				for Nms_a in range(0, self.Nms_tot(Nmm)):
					ave_Nmm += Nmm * self.boltzmann(Nmm, Nmm_a, Nms_a)
					Z       += self.boltzmann(Nmm, Nmm_a, Nms_a)

		return (ave_Nmm, Z)

	def get_Nmm(self):

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
		
		# make the pool
		pool = multiprocessing.Pool(processes=self.nproc)
		
		# run the parallel calc
		results = pool.starmap(self.run_Nmm_sum, zip(segments)) 

		ave_Nmm = 0
		Z       = 0

		for res in results:
			ave_Nmm += res[0]
			Z       += res[1]

		return ave_Nmm/Z

	def plot_Nmm(self):
		aNmm_store = list()

		for T in self.T_range:
			print(f"@T = {T}...", flush=True)
			self.T = T 
			ave_Nmm = self.get_Nmm()
			aNmm_store.append(ave_Nmm)

		# set up the plot 
		fig = plt.figure(figsize=(2, 2))
		ax  = plt.axes()

		# We change the fontsize of minor ticks label 
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.plot(T_range, aNmm_store, marker='o', mec='k', c='steelblue', ls='--', lw=2, markersize=6)
		ax.minorticks_on()
		fig.savefig(args.o, dpi=1200, bbox_inches="tight")

		f = open(args.csv, 'w')
		f.write(f"Temperature | Nmm\n")
		for i in range(len(self.T_range)):
			f.write(f"{self.T_range[i]} | {aNmm_store[i]}\n")
		f.close()

		return 

	# this is to get Nms aligned 
	def run_Nmsa_sum(self, Nmm_range):
		ave_Nmsa = 0
		Z        = 0

		for Nmm in Nmm_range:
			for Nmm_a in range(0, Nmm):
				for Nms_a in range(0, self.Nms_tot(Nmm)):
					ave_Nmsa += Nms_a * self.boltzmann(Nmm, Nmm_a, Nms_a)
					Z        += self.boltzmann(Nmm, Nmm_a, Nms_a)
		
		return (ave_Nmsa, Z)

	def get_Nmsa(self):

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
		
		# make the pool
		pool = multiprocessing.Pool(processes=self.nproc)
		
		# run the parallel calc
		results = pool.starmap(self.run_Nmsa_sum, zip(segments)) 

		ave_Nmsa = 0
		Z       = 0

		for res in results:
			ave_Nmsa += res[0]
			Z        += res[1]

		return ave_Nmsa/Z

	def plot_Nmsa(self):
		aNmsa_store = list()

		for T in self.T_range:
			print(f"@T = {T}...", flush=True)
			self.T = T 
			ave_Nmsa = self.get_Nmsa()
			aNmsa_store.append(ave_Nmsa)

		# set up the plot 
		fig = plt.figure(figsize=(2, 2))
		ax  = plt.axes()

		# We change the fontsize of minor ticks label 
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.tick_params(axis='both', which='major', labelsize=6)
		ax.plot(T_range, aNmsa_store, marker='o', mec='k', c='springgreen', ls='--', lw=2, markersize=6)
		ax.minorticks_on()
		fig.savefig(args.o, dpi=1200, bbox_inches="tight")

		f = open(args.csv, 'w')
		f.write(f"Temperature | Nmsa\n")
		for i in range(len(self.T_range)):
			f.write(f"{self.T_range[i]} | {aNmsa_store[i]}\n")
		f.close()

		return 

if __name__=="__main__":

	# get some time counters 
	start_global = time.time()

	# define properties of objects
	Nm      = 32
	T_range = args.T
	nproc   = args.nproc
	maxwell = Maxwell(Nm, T_range, nproc)

	if args.style == "distribution":
		maxwell.plot_distributions()

	elif args.style == "energy":
		maxwell.plot_energy()

	elif args.style == "free_energy":
		maxwell.plot_free_energy()

	elif args.style == "mm_contacts":
		maxwell.plot_Nmm()

	elif args.style == "aligned_ms_contacts":
		maxwell.plot_Nmsa()

	elif args.style == "heat_capacity":
		maxwell.plot_variance()

	else:
		print(f"You did not enter a style. Exiting...", flush=True)
		exit()

	# get the timings
	stop_global = time.time()
	print(f"Time for computation is {stop_global-start_global} seconds.", flush=True)
