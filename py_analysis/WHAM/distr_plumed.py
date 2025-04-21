#!/home/satyend/.conda/envs/phase/bin/python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pandas as pd
import re
import os
import argparse

parser = argparse.ArgumentParser(description="Plot distributions from COLVAR files.")
parser.add_argument ("-i", dest='i', type=str,   action='store', help="Name of image to make.")
args = parser.parse_args ()

if __name__=="__main__":

	fig    = plt.figure(num=0)
	ax     = plt.axes()
	n_bins = 100

	norm     = plt.Normalize(vmin=9, vmax=20)
	colormap = mpl.colormaps['viridis']

	files = [f for f in os.listdir(".") if f.startswith("AT_")]
	files = sorted(files, key=lambda x: float(x.split('_')[1]))
	cleaned_files = list()

	for fil in files:
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR")]
		if len(colvar_files) != 0:
			cleaned_files.append(fil)
	
	print("Go through all the files.")
	for fil in cleaned_files:
		# print(f"In {fil}.", flush=True)
		rg = np.array([])
		c = colormap(norm(float(fil[3:])))
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR")]
		for idx,cf in enumerate(colvar_files):
			df = pd.read_csv(fil+"/"+cf, skiprows=1, engine='python', names=["time", "rg", "bias"], sep='\s+')
			if idx == 0:
				rg = np.hstack((rg, df["rg"].values[-700:]))
			else:
				rg = np.hstack((rg, df["rg"].values))
		counts, bins = np.histogram(rg, bins=n_bins, range=[9, 18])
		bin_centers  = 0.5 * (bins[1:] + bins[:-1])
		ax.hist(rg, bin_centers, range=[9, 18], alpha=1, label=f"$x_0 = {fil[3:]}$", color='white', edgecolor=c)

	ax.legend()
	ax.set_xlim(8, 19)
	ax.set_xticks(np.arange(8,20))
	fig.savefig(args.i+".png", dpi=1200, bbox_inches="tight")

	fig    = plt.figure(num=1)
	ax     = plt.axes()

	print("Go through some of the files.")
	for fil in cleaned_files[::2]:
		rg = np.array([])
		c = colormap(norm(float(fil[3:])))
		colvar_files = [f for f in os.listdir(fil) if f.startswith("COLVAR")]
		for cf in colvar_files:
			df = pd.read_csv(fil+"/"+cf, skiprows=1, engine='python', names=["time", "rg", "bias"], sep='\s+')
			if idx == 0:
				rg = np.hstack((rg, df["rg"].values[-700:]))
			else:
				rg = np.hstack((rg, df["rg"].values))
		counts, bins = np.histogram(rg, bins=n_bins, range=[9, 18])
		bin_centers  = 0.5 * (bins[1:] + bins[:-1])
		# ax.plot(bin_centers, counts, color=c, lw=1, label=f"$x_0 = {fil[3:]}$")
		ax.hist(rg, bin_centers, range=[9, 18], alpha=1, label=f"$x_0 = {fil[3:]}$", color='white', edgecolor=c)

	ax.legend()
	ax.set_xlim(8, 19)
	ax.set_xticks(np.arange(8,20))
	fig.savefig(args.i+"_cut.png", dpi=1200, bbox_inches="tight")
