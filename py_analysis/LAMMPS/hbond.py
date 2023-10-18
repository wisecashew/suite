#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import argparse
import time

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the hydrogen bonds over time.")
parser.add_argument('--coords',    dest='c',     action='store', type=str,   help='Name of the coordinate dump.', default=None)
parser.add_argument('--data',      dest='data',  action='store', type=str,   help='Name of the data file.', default=None)
parser.add_argument('--image',     dest='img',   action='store', type=str,   help='Name of image to be created.', default=None)
parser.add_argument('--hydrogens', dest='h',     action='store', type=str,   help='Identification for hydrogens.', default=None)
parser.add_argument('--acceptors', dest='acc',   action='store', type=str,   help='Identification for acceptor atoms.', default=None)
parser.add_argument('--donors',    dest='don',   action='store', type=str,   help='Identification for donor atoms.', default=None)
args = parser.parse_args()


if __name__=="__main__":

	start = time.time()

	try:
		u = mda.Universe(args.data, args.c, format="LAMMPS", lengthunit="angstrom", timeunit="ns")
	except:
		u = mda.Universe(args.data, args.c, atom_style="id type x y z", lengthunit="angstrom", timeunit="fs") #, lengthunit="angstrom", timeunit="ns")

	# access atom types for hbonding
	hbonds = HBA(universe=u, hydrogens_sel=args.h, acceptors_sel=args.acc, donors_sel=args.don)
	hbonds.run()

	times  = hbonds.times
	counts = hbonds.count_by_time()

	fig = plt.figure(figsize=(2,2))
	ax  = plt.axes()

	ax.plot(times/1e+6, counts, lw=1, ls='--', color='steelblue', marker='o', markersize=2, mec='k')
	fig.savefig(args.img, bbox_inches="tight", dpi=1200)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)
