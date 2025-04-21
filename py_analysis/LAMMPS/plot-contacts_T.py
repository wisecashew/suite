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

	df = pd.read_csv(args.df, engine="python", names=['T', "Nmm_mean", "Nmm_err", "Nms_mean", "Nms_err", "Nmsa_mean", "Nmsa_err", "Ns_mean", "Ns_err"], sep='\s+')

	# figure 0
	cols = ["coral", "steelblue", "lavender", "forestgreen"]
	fig  = plt.figure(num=0, figsize=(3,3))
	ax   = plt.axes() 
	ax.errorbar(df['T'].values, df["Nmm_mean"].values, yerr=df["Nmm_err"].values, c=cols[0], marker='o', ls='--', lw=1, markersize=4, mec='k')
	ax.set_xlim(np.min(df['T'].values), np.max(df['T'].values))
	ax.set_ylim(70, 90)
	ax.axhline(y=80)
	fig.savefig(args.img+"_Nmm", dpi=1200, bbox_inches="tight")
	plt.close(fig)

	# figure 1
	fig  = plt.figure(num=1, figsize=(3,3))
	ax   = plt.axes() 
	ax.errorbar(df['T'].values, df["Nms_mean"].values, yerr=df["Nms_err"].values, c=cols[1], marker='o', ls='--', lw=1, markersize=4, mec='k')
	ax.set_xlim(np.min(df['T'].values), np.max(df['T'].values))
	fig.savefig(args.img+"_Nms", dpi=1200, bbox_inches="tight")
	plt.close(fig)

	# figure 2
	fig  = plt.figure(num=2, figsize=(3,3))
	ax   = plt.axes() 
	ax.errorbar(df['T'].values, df["Nmsa_mean"].values, yerr=df["Nmsa_err"].values, c=cols[1], marker='o', ls='--', lw=1, markersize=4, mec='k')
	ax.set_xlim(np.min(df['T'].values), np.max(df['T'].values))
	fig.savefig(args.img+"_Nmsa", dpi=1200, bbox_inches="tight")
	plt.close(fig)

	# figure 3
	fig  = plt.figure(num=2, figsize=(3,3))
	ax   = plt.axes() 
	ax.errorbar(df['T'].values, df["Ns_mean"].values, yerr=df["Ns_err"].values, c=cols[1], marker='o', ls='--', lw=1, markersize=4, mec='k')
	ax.set_xlim(np.min(df['T'].values), np.max(df['T'].values))
	fig.savefig(args.img+"_Ns", dpi=1200, bbox_inches="tight")
	plt.close(fig)
