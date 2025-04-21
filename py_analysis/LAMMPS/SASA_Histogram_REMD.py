#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the SASA histogram.")
parser.add_argument("--image",     dest="img",    action="store", type=str,   default=None, required=True,  help="Name of image to be created.")
parser.add_argument("--nchunk",    dest="nchunk", action="store", type=int,   default=1,    required=False, help="Number of chunks to make.")
parser.add_argument("--avg",       dest="avg",    action="store", type=str,   default=None, required=True,  help="Enter name of file to dump averages to.")
parser.add_argument("--yulim",     dest="yulim",  action="store", type=float, default=None, required=False, help="Enter upper limit along y-axis.")
args = parser.parse_args()

def chunked_histogram(data, num_chunks, bins, hrange):

	# split the data into chunks
	chunk_size = len(data) // num_chunks
	print(f"chunk_size = {chunk_size}")
	print(f"num_chunks = {num_chunks}")
	chunks = [data[i * chunk_size:(i+1) * chunk_size] for i in range(num_chunks)]

	# initialize an empty list to hold histograms
	histograms = []

	# compute histogram for each chunk
	for chunk in chunks:
		hist, bin_edges = np.histogram(chunk, bins=bins, range=hrange)
		histograms.append(hist)

	# convert the lsit of histograms into a numpy array for easier manipulation
	histograms = np.array(histograms)

	# calculate the mean and standard deviation for each bin
	mean_hist  = np.mean(histograms, axis=0)
	error_hist = np.std(histograms, axis=0)/np.sqrt(num_chunks)

	# calculate the bin centers for plotting
	bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

	return bin_centers, mean_hist, error_hist

def block_averaging(data, num_chunks):

	# split the data
	chunk_size = len(data) // num_chunks 
	chunks     = [data[i * chunk_size:(i+1) * chunk_size] for i in range(num_chunks)]
	means      = [np.mean(chunk) for chunk in chunks]

	average    = np.mean(means)
	error      = np.std(means)/num_chunks

	return average, error

if __name__=="__main__":

	# track how the temperature index goes around
	# indices go from 0 to 39
	T_range   = np.linspace(275, 360, 40)
	SASA_over_T = dict()
	for idx in range(0, 40):
		SASA_over_T[idx] = []

	DF_list = []
	for idx in range(0, 40):
		df = pd.read_csv("TOT_SASA_FILES/sasa."+str(idx)+".dump", sep='\s+', names=["timestep", "sasa"], engine="python", skiprows=200)
		DF_list.append(df)
	
	print(f"Temperature = ", end='', flush=True)
	for Tidx in range(0, 40):
		if Tidx < 39:
			print(f"{T_range[Tidx]:.2f}, ",end='', flush=True)
		else:
			print(f"{T_range[Tidx]:.2f}.",end='', flush=True)
		output_file = open("TOT_OUTPUT_FILES/output.net", 'r')
		for line in output_file:
			if re.match(r"^[0-9]+", line):
				line       = line.strip().split()
				line       = [int(i) for i in line]
				ts         = line[0]
				indices    = line[1:]
				replica_id = indices.index(Tidx)
				timesteps  = np.array(DF_list[replica_id]["timestep"].values, dtype=float)
				if np.sum(timesteps==ts) == 1:
					SASA_over_T[Tidx].extend(DF_list[replica_id]["sasa"].values[timesteps==ts].tolist())
		output_file.close()

	print()
	print(f"Make images...", flush=True)

	num_chunks = args.nchunk
	avg_dumpfile = open(args.avg, 'w')
	hrange     = [3000, 7000]
	for Tidx in range(0, 40):
		print(f"@ T = {T_range[Tidx]:.2f} and ", end='')
		fig = plt.figure(num=Tidx, figsize=(3,3))
		ax  = plt.axes()
		print(f"number of samples = {len(SASA_over_T[Tidx])}.")
		bin_centers, mean_hist, error_hist = chunked_histogram(SASA_over_T[Tidx], num_chunks, 50, hrange)
		ax.errorbar(bin_centers, mean_hist, yerr=error_hist, marker='o', markeredgecolor='k', c="coral")
		mean_sasa, err_sasa = block_averaging(SASA_over_T[Tidx], num_chunks)
		avg_dumpfile.write(f"{T_range[Tidx]:.2f} {mean_sasa:.2f} {err_sasa:.2f}\n")
		ax.set_title(f"T = {T_range[Tidx]:.2f} K")
		ax.set_xlabel("SASA (angstrom^2)")
		ax.set_ylabel("Frequency")
		ax.set_xlim(2500, 7500)
		yulim = 10**(np.ceil(np.log10(np.max(mean_hist)))) if args.yulim is None else args.yulim
		ax.set_ylim(0, yulim)
		ax.set_xticks      ([2500, 3500, 4500, 5500, 6500, 7500])
		ax.set_xticklabels ([2500, 3500, 4500, 5500, 6500, 7500])
		fig.savefig("TOT_IMAGES/SASA/"+args.img+f"_n{num_chunks}_{T_range[Tidx]:.2f}.png", dpi=1200, bbox_inches="tight")
		plt.close(fig)
	avg_dumpfile.close()
