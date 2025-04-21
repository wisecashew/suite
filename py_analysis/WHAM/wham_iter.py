#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import re
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Run the WHAM algorithm.")
parser.add_argument ("-i", dest='i', type=str,   action='store', help="Name of image to make.")
parser.add_argument ("-T", dest='T', type=float, action='store', help="Temperature of simulation.")
parser.add_argument("--half", dest='h', action='store_true', default=False, help="Only consider half the data.")
args = parser.parse_args ()

# get the boltzmann constant
kb = 0.0019872041
T  = args.T
beta = 1/ (kb * T)
print(f"Thermodynamic beta = {beta}")

# get the energetic parameters
def read_plumed(pl_file):

	# define the regex patterns
	pattern_kappa = r"KAPPA=([\d\.]+)"
	pattern_at    = r"AT=([\d\.]+)"

	# instantiate 
	kappa = None
	at    = None 

	# go through the files
	with open(pl_file, 'r') as f:
		for line in f:
			if kappa is None:
				match = re.findall(pattern_kappa, line)
				if match:
					kappa = float(match[0])
			if at is None:
				match = re.findall(pattern_at, line)
				if match:
					at = float(match[0])
			if (not (kappa is None)) and (not (at is None)):
				break
	return kappa, at

# set up the histograms and stuff
def extract_info(cleaned_files, n_bins):
	H         = list()
	N_samples = list()
	U_B       = list()

	# set up the Rg container
	rg = np.array([])

	for fil in cleaned_files:
		kappa, x0 = read_plumed(fil + "/plumed.dat")
		U_B.append(lambda Rg, k=kappa, x=x0: 1/2 * k * (Rg - x) ** 2)

		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR") and f.endswith("run")]
		colvar_files = sorted(colvar_files, key=lambda x:float(x.split('.')[1]))
		for cf in colvar_files:
			df = pd.read_csv(fil+"/"+cf, skiprows=1, engine="python", names=["time", "rg", "bias"], sep='\s+')
			rg = np.hstack((rg, df["rg"].values))
	
	# get the minima and maxima
	rg_min = np.min(rg)
	rg_max = np.max(rg)

	# print(f"rg_min = {rg_min}, rg_max = {rg_max}", flush=True)

	# get the bins
	bins        = np.linspace(rg_min, rg_max, n_bins+1)
	bin_centers = (bins[1:] + bins[:-1])/2

	# run through the cleaned files
	for fil in cleaned_files:
		# set up the Rg container
		rg = np.array([])

		# get the colvar files
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR") and f.endswith("run")]
		for idx,cf in enumerate(colvar_files):
			df = pd.read_csv(fil+"/"+cf, skiprows=1, engine="python", names=["time", "rg", "bias"], sep='\s+')
			if idx == 0:
				rg = np.hstack((rg, df["rg"].values[-700:]))
			else:
				rg = np.hstack((rg, df["rg"].values[10:]))

		counts = np.histogram(rg, bins=bins, range=[rg_min, rg_max])[0]
		H.append(counts)
		N_samples.append(len(rg))

	return H, N_samples, U_B, bins, bin_centers

if __name__=="__main__":

	# set up the WHAM binning parameters
	n_bins = 100

	# grab the data
	sim_files     = [f for f in os.listdir(".") if f.startswith("AT_")]
	sims_files    = sorted(sim_files, key=lambda x: float(x.split('_')[1]))
	cleaned_files = list()

	# get the files with the relevant information
	for fil in sims_files:
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR")]
		if len(colvar_files) != 0:
			cleaned_files.append(fil)

	# get all the info
	histograms, N_samples, U_B, bins, bin_centers = extract_info(cleaned_files, n_bins)

	# get some data 
	n_sims = len(U_B)
	print(f"Number of sims is {n_sims}.", flush=True)
	Z      = np.ones(n_sims) # np.random.uniform(low=1, high=10, size=(n_sims,))
	Z_new  = np.ones(n_sims)
	pj     = np.zeros_like(bin_centers)

	# set up some convergence criteron
	tol     = 1e-6
	maxiter = 1000
	converged = False

	if args.h:
		jump = 2
	else:
		jump = 1

	it = 0
	while not converged and it < maxiter:
		# set up the iterative solve
		for ibc, bc in enumerate(bin_centers):
			numer   = np.sum( [histograms[s][ibc] for s in range(0, n_sims, jump)] )
			denom   = np.sum( [N_samples[s]*Z[s]*np.exp(-beta * U_B[s](bc)) for s in range(0, n_sims, jump)] )
			pj[ibc] = numer/denom

		# update the free energies 
		for s in range(0, n_sims, jump):
			prob_mass = 0
			for ibc, bc in enumerate(bin_centers):
				prob_mass += np.exp(-beta * U_B[s](bc)) * pj[ibc]
			Z_new[s] = 1/prob_mass
		
		# get the error, update iterator
		diff = np.linalg.norm(Z-Z_new)
		# print(f"@iteration {it}: error = {diff}", flush=True)
		if diff > tol:
			Z = Z_new
		else:
			converged = True
		it += 1

	# The free energy surface is related to -log(weight)
	free_energy_sim     = -kb * T * np.log(Z)
	free_energy_sim    -= np.min(free_energy_sim)
	free_energy_window  = -kb * T * np.log(pj)
	free_energy_window -= np.min(free_energy_window)

	# get the figure set up
	fig = plt.figure(num=0, figsize=(3,3))
	ax  = plt.axes()

	# plot the WHAM'd surface
	ax.plot(bin_centers, free_energy_window, lw=1, marker='o', mec='k', label="WHAM")
	# print(f"Minima @ {bin_centers[np.argmin(free_energy_window)]}")
	ax.legend()
	ax.set_xlim(9, 18)
	ax.set_ylim(0, 6)
	ax.set_xticks(np.arange(9, 19, 1))
	ax.set_xticklabels(np.arange(9, 19, 1))
	ax.minorticks_on()

	# save the figure
	fig.savefig(args.i, dpi=1200, bbox_inches="tight")
