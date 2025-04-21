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
args = parser.parse_args ()

# get the boltzmann constant
kb = 0.0019872041

# do WHAM
def wham(bias,
		*,
		frame_weight=None,
		traj_weight=None,
		T: float = 1.0,
		maxiter: int = 1000,
		threshold: float = 1e-40,
		verbose: bool = False):

	nframes = bias.shape[0]
	ntraj   = bias.shape[1]

	# default values
	if frame_weight is None:
		frame_weight = np.ones(nframes)
	if traj_weight is None:
		traj_weight = np.ones(ntraj)

	assert len(traj_weight)  == ntraj
	assert len(frame_weight) == nframes

	# divide by T once for all
	shifted_bias = bias/T

	# track shifts
	shifts0       = np.min(shifted_bias, axis=0)
	shifted_bias -= shifts0[np.newaxis,:]
	shifts1       = np.min(shifted_bias, axis=1)
	shifted_bias -= shifts1[:,np.newaxis]

	# do exponentials only once
	expv = np.exp(-shifted_bias)

	# instantiate the partition function
	Z = np.ones(ntraj)

	# copy it over
	Zold = Z

	if verbose:
		sys.stderr.write("WHAM: start\n")
	for nit in range(maxiter):
		# find unnormalized weights
		weight = 1.0/np.matmul(expv, traj_weight/Z)*frame_weight
		# update partition functions
		Z = np.matmul(weight, expv)
		# normalize the partition functions
		Z /= np.sum(Z*traj_weight)
		# monitor change in partition functions
		eps = np.sum(np.log(Z/Zold)**2)
		Zold = Z
		if verbose:
			sys.stderr.write("WHAM: iteration "+str(nit)+" eps "+str(eps)+"\n")
		if eps < threshold:
			break
	nfev=nit
	logW = np.log(weight) + shifts1

	if verbose:
		sys.stderr.write("WHAM: end")

	return {"logW":logW, "logZ":np.log(Z)-shifts0, "nit":nit, "eps":eps}

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

	print(f"cleaned_files = {cleaned_files}")
	for fil in cleaned_files:
		kappa, x0 = read_plumed(fil + "/plumed.dat")
		Usim      = lambda Rg, k=kappa, x=x0: 1/2 * k * (Rg - x) ** 2
		U_B.append(Usim)

		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR") and f.endswith("run")]
		colvar_files = sorted(colvar_files, key=lambda x:float(x.split('.')[1]))
		print(colvar_files)
		for cf in colvar_files:
			df = pd.read_csv(fil+"/"+cf, skiprows=1, engine="python", names=["time", "rg", "bias"], sep='\s+')
			rg = np.hstack((rg, df["rg"].values))
	
	# get the minima and maxima
	rg_min = np.min(rg) 
	rg_max = np.max(rg)

	# get the bins
	bins        = np.linspace(rg_min, rg_max, n_bins+1)
	bin_centers = (bins[1:] + bins[:-1])/2

	for fil in cleaned_files:
		kappa, x0 = read_plumed(fil + "/plumed.dat")
		Usim      = lambda Rg: 1/2 * kappa * (Rg - x0) ** 2
		U_B.append(Usim)

		# set up the Rg container
		rg = np.array([])

		# get the colvar files
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR")]
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
	F_i    = np.zeros(n_sims)

	# set up the variables to run the solver
	bias = np.zeros((n_bins, n_sims))

	print(f"Getting the matrix.")
	for i, x in enumerate(bin_centers):
		for ns in range(n_sims):
			bias[i, ns] = U_B[ns](x)

	# get the histograms from biased simulations
	frame_weight = np.sum(histograms, axis=0)
	results = wham(bias, frame_weight=frame_weight, T=kb*args.T, maxiter=1000, verbose=True)

	# Extract the reconstructed free energy surface
	logW = results['logW']

	# The free energy surface is related to -log(weight)
	free_energy_surface  = -logW
	free_energy_surface -= np.min(free_energy_surface)

	# get the figure set up
	fig = plt.figure(num=0, figsize=(3,3))
	ax  = plt.axes()

	# plot the WHAM'd surface
	ax.plot(bin_centers, free_energy_surface, lw=1, marker='o', mec='k', label="WHAM")
	print(f"Minima @ {bin_centers[np.argmin(free_energy_surface)]}")
	ax.legend()

	# save the figure
	fig.savefig(args.i, dpi=1200, bbox_inches="tight")
