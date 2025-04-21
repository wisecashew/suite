#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the Rg histogram.")
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
	contacts_over_T = dict()
	for idx in range(0, 40):
		contacts_over_T[idx] = [] # list()

	DF_list = [] # list()
	print(f"Going to compile \"contacts.n.dump\".", flush=True)
	for idx in range(0, 40):
		df = pd.read_csv("TOT_DRY_LATT/contacts."+str(idx)+".dump", sep='\s+', names=["timestep", "Nmm", "Nms"], engine="python", skiprows=2000)
		DF_list.append(df)
	
	print(f"Temperature = ", end='', flush=True)
	for Tidx in range(0, 40):
		if Tidx < 39:
			print(f"{T_range[Tidx]:.2f}, ", end='', flush=True)
		else:
			print(f"{T_range[Tidx]:.2f}." , end='', flush=True)
		output_file = open("TOT_OUTPUT_FILES/output.net", 'r')
		for line in output_file:
			if re.match(r"^[0-9]+", line):
				line       = line.strip().split()
				line       = [int(i) for i in line]
				ts         = line[0]
				indices    = line[1:]
				replica_id = indices.index(Tidx)
				timesteps  = np.array(DF_list[replica_id]["timestep"].values, dtype=int)
				if np.sum(timesteps==ts) == 1:
					contacts_over_T[Tidx].extend(DF_list[replica_id]["Nmm"].values[timesteps==ts].tolist())
		output_file.close()

	print()
	print(f"Make images...", flush=True)


	num_chunks = args.nchunk
	avg_dumpfile = open(args.avg, 'w')
	hrange     = [30, 450]
	for Tidx in range(0, 40):
		print(f"@ T = {T_range[Tidx]:.2f} and ", end='')
		fig = plt.figure(num=Tidx, figsize=(3,3))
		ax  = plt.axes()
		print(f"number of samples = {len(contacts_over_T[Tidx])}.")
		bin_centers, mean_hist, error_hist = chunked_histogram(contacts_over_T[Tidx], num_chunks, 50, hrange)
		mean_contacts, err_contacts = block_averaging(contacts_over_T[Tidx], num_chunks)
		avg_dumpfile.write(f"{T_range[Tidx]:.2f} {mean_contacts:.2f} {err_contacts:.2f}\n")
		ax.errorbar(bin_centers, mean_hist, yerr=error_hist, marker='o', markeredgecolor='k', c="steelblue")
		ax.axvline(x=80, ls="--", c="darkred")
		ax.set_title(f"T = {T_range[Tidx]:.2f} K")
		ax.set_xlabel("$N_{mm}$")
		ax.set_ylabel("Frequency")
		# ax.set_xlim(5, 30)
		# yulim = 10**(np.ceil(np.log10(np.max(mean_hist)))) if args.yulim is None else args.yulim
		# ax.set_ylim(0, yulim)
		# ax.set_xticks([5, 10, 12.5, 15, 20, 25, 30])
		# ax.set_xticklabels([5, 10, 12.5, 15, 20, 25, 30])
		fig.savefig("TOT_IMAGES/NMM/"+args.img+f"n-{num_chunks}_{T_range[Tidx]:.2f}.png", dpi=1200, bbox_inches="tight")
		plt.close(fig)
	avg_dumpfile.close()
