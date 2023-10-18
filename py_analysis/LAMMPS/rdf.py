#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the flory exponent from that file.")
parser.add_argument('--coords',    dest='c',     action='store', type=str,   help='Name of energy dump file to parse information.', default=None)
parser.add_argument('--image',     dest='img',   action='store', type=str,   help='Name of image to be created.', default=None)
parser.add_argument('--atm1',      dest='atm1',  action='store', type=str,   help='First atom type to be involved in rdf calculation.' , default=None)
parser.add_argument('--atm2',      dest='atm2',  action='store', type=str,   help='Second atom type to be involved in rdf calculation.', default=None)
parser.add_argument('--rdf-range', dest='range', action='store', type=float, help='Enter range of RDF calculation.', default=None)
args = parser.parse_args()


if __name__=="__main__":

	coords = mda.coordinates.LAMMPS.DumpReader(args.c)
	u = mda.Universe(args.c, atom_style="id type x y z")

	# access atom types
	atm1 = u.select_atoms(args.atm1)
	atm2 = u.select_atoms(args.atm2)
	if args.atm1 == args.atm2:
		irdf = rdf.InterRDF(atm1, atm2, nbins=100, range=(0.0, args.range), exclusion_block=(1,1))
	else:
		irdf = rdf.InterRDF(atm1, atm2, nbins=100, range=(0.0, args.range))
	irdf.run()

	plt.plot(irdf.bins, irdf.rdf, c='steelblue')
	if args.img == None:
		plt.savefig("rdf", dpi=1200, bbox_inches="tight")
	else:
		plt.savefig(args.img, dpi=1200, bbox_inches="tight")


