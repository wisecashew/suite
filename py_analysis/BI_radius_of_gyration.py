#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Analyze radius of gyration.")
parser.add_argument("-i", dest='i', type=str, action='store', help="Name of rg file.")
args = parser.parse_args()


if __name__=="__main__":

	fig = plt.figure ()
	ax  = plt.axes   ()

	df = pd.read_csv (args.i, engine='python', sep='\s', names=["time", "Rg"], skiprows=2)
	Rg = df["Rg"].values
	(n_h, bins_h, patches) = ax.hist (Rg, density=True, bins=100)
	fig.savefig ("rg_hist", dpi=1200, bbox_inches="tight")

	bin_collect = (bins_h[0:-1]+bins_h[1:])/2
	F_x1 = -2/3 * np.log (n_h)

	fig1 = plt.figure ()
	ax1  = plt.axes   ()

	ax1.plot (bin_collect, F_x1, label="unbiased")
	fig1.savefig ("rg_fes", dpi=1200, bbox_inches="tight")
