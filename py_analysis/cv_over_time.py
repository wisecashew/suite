#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import re

parser = argparse.ArgumentParser(description="Plot the FES from COLVARS.")
parser.add_argument("-i", dest='i', type=str, action='store', help="Enter address of CV trajectory.")
parser.add_argument("--xllim", dest='xllim', type=float, action='store', help="Enter lower x limit.", default=None)
parser.add_argument("--xulim", dest='xulim', type=float, action='store', help="Enter upper x limit.", default=None)
parser.add_argument("--xlogspace", dest='xls', action='store_true', help="Enter if you want x-axis logspaced.", default=False)
parser.add_argument("--yllim", dest='yllim', type=float, action='store', help="Enter lower y limit.", default=None)
parser.add_argument("--yulim", dest='yulim', type=float, action='store', help="Enter upper y limit.", default=None)
parser.add_argument("--img", dest='o', type=str, action='store', help="Enter address of image to be made.")


args   = parser.parse_args()

if __name__=="__main__":

	# print (list_of_files)
	fig = plt.figure(figsize=(3,3))
	ax  = plt.axes()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')

	df = pd.read_csv(args.i, sep='\s+', names=["time", "rg", "restraint"], engine='python', comment='#')
	time = list(range(len(df["time"].values)))
	# print(time)
	ax.scatter(time, df["rg"].values, c="steelblue", s=0.5)

	if args.yllim == None:
		ax.set_ylim(bottom=np.min(df["rg"].values)*0.9)
	else:
		ax.set_ylim(bottom=args.yllim)

	if args.yllim == None:
		ax.set_ylim(top=np.max(df["rg"].values)*1.1)
	else:
		ax.set_ylim(top=args.yulim)

	if args.xllim == None:
		ax.set_xlim(left=time[0])
	else:
		ax.set_xlim(left=args.xllim)


	if args.xulim == None:
		ax.set_xlim(right=time[-1])
	else:
		ax.set_xlim(right=args.xulim)

	if args.xls:
		if args.xllim == None:
			ax.set_xlim(left=1)
		else:
			ax.set_xlim(left=args.xllim)
		ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[-1], 6))
		ax.set_xticklabels(['']*len(ax.get_xticklabels()))
		ax.set_xticklabels(np.linspace(ax.get_xlim()[0], ax.get_xlim()[-1], 6, dtype=int))

	else:
		ax.minorticks_on()
		ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

	ax.set_xlabel("time (ps)")
	ax.set_ylabel("$R_g$ (A)")

	if (".png" in args.o[-4:]):
		img_name = args.o
	elif ("." in args.o):
		img_name = args.o + ".png"
	else:
		img_name = args.o
	fig.savefig (img_name, dpi=1200, bbox_inches="tight")

