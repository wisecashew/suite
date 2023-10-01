#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import re

parser = argparse.ArgumentParser(description="Plot the FES from `plumed sum_hills`.")
parser.add_argument("-i", dest='i', type=str, action='store', help="Enter address of FES file.")
parser.add_argument("--xllim", dest='xllim', type=float, action='store', help="Enter lower x limit.", default=None)
parser.add_argument("--xulim", dest='xulim', type=float, action='store', help="Enter upper x limit.", default=None)
parser.add_argument("--yllim", dest='yllim', type=float, action='store', help="Enter lower y limit.", default=None)
parser.add_argument("--yulim", dest='yulim', type=float, action='store', help="Enter upper y limit.", default=None)
parser.add_argument("--img", dest='o', type=str, action='store', help="Enter address of image to be made.")
parser.add_argument("--multiple", dest='multiple', action='store_true', help="Enter name of file over which the evolution is taking place.", default=False)
parser.add_argument("--address", dest='addr', action='store', type=str, help="Enter address of fes evolution.", default=".")
parser.add_argument("--start", dest='start', action='store', type=int, help="Enter index at which you need to start probing.", default=0)
parser.add_argument("--stop", dest='stop', action='store', type=int, help="Enter index at which you need to stop probing.")
parser.add_argument("--skip", dest='skip', action='store', type=int, help="Enter how many files to skip (this option only matters if --multiple is called).")

def get_fes_files(address):
	"""Returns a list of all the files in the current directory with names "fes_xxx.dat" where "xxx" is a number, sorted by the number "xxx".

	Returns:
	A list of strings, where each string is the name of a file.
	"""

	fes_files = []
	for filename in os.listdir(address):
		if filename.startswith("fes_") and filename.endswith(".dat"):
			try:
				int(filename[4:-4])
				fes_files.append(filename)
			except ValueError:
				pass

	# Sort the files by the number "xxx".
	fes_files.sort(key=lambda filename: int(filename[4:-4]))

	return fes_files


args   = parser.parse_args()

if __name__=="__main__":

	if args.multiple:
		list_of_files = get_fes_files(args.addr)
		# print (list_of_files)
		fig = plt.figure(figsize=(3,3))
		ax  = plt.axes()
		ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
		for idx, f in enumerate(list_of_files[args.start:args.stop:args.skip]):
			df = pd.read_csv (args.addr+f, engine="python", skiprows=5, names=["rg", "fes", "dfes"], sep='\s+')
			ax.plot (df["rg"].values, df["fes"].values, lw=1, label=f"iter ={args.start+idx*args.skip}")

		ax.legend (loc="upper right", fontsize=2)

	else:
		df = pd.read_csv (args.i, engine="python", skiprows=5, names=["rg", "fes", "dfes"], sep='\s+')
		fig = plt.figure(figsize=(3,3))
		ax  = plt.axes ()
		ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
		ax.plot (df["rg"].values, df["fes"].values, c='steelblue', lw=1)

	if args.yllim == None or args.yulim == None:
		ax.set_ylim((0, 250))
	else:
		ax.set_ylim((args.yllim, args.yulim))

	if args.xllim == None or args.xulim == None:
		ax.set_xlim((0, 30))
	else:
		ax.set_xlim((args.xllim, args.xulim))


	if "." in args.o:
		img_name = args.o + ".png"
	else:
		img_name = args.o
	ax.minorticks_on()
	fig.savefig (img_name, dpi=1200, bbox_inches="tight")
