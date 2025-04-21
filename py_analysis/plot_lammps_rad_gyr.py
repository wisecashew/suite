#!/home/satyend/.conda/envs/phase/bin/python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description = "Read a LAMMPS .avg file for gyration.")
parser.add_argument('-i', metavar='rad_gyr.avg', dest='i', action='store', help='enter address of xvg file', type=str)
args = parser.parse_args()

if __name__=="__main__":

	f = open(args.i, 'r')
	df = pd.read_csv(f, names=["t", "Rg", "Rgx", "Rgy", "Rgz"], delim_whitespace=True, comment='#')
	f.close()

	print ("Done reading file.", flush=True)
	print ("Plotting.", flush=True)

	rg_vals = df["Rg"].values
	plt.scatter(range(len(df["t"])), rg_vals, c='steelblue', s=4)
	plt.savefig ("rg_dist", dpi=1200)

	print("Mean Rg is {:0.2f}".format(np.mean(rg_vals)))
	print("Standard deviation of Rg is {:0.2f}".format(np.std(rg_vals)))

