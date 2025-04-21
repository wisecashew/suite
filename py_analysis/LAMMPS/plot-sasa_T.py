#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot sasa over temperature from a simple 2-column file.")
parser.add_argument("--dumpfile", dest="df",  type=str, action="store", required=True, help="enter address of dump file.")
parser.add_argument("--img",      dest="img", type=str, action="store", required=True, help="enter name of image to create.")
args = parser.parse_args()

if __name__=="__main__":

	df = pd.read_csv(args.df, engine="python", names=['T', "sasa", "error"], sep='\s+')
	print(df)
	fig, ax = plt.figure(num=0, figsize=(3,3)), plt.axes()
	
	ax.errorbar(df['T'].values, df["sasa"].values, yerr=df["error"].values, c="coral", marker='o', ls='--', lw=1, markersize=4, mec='k')
	ax.set_xlim(np.min(df['T'].values), np.max(df['T'].values))

	fig.savefig(args.img, dpi=1200, bbox_inches="tight")

