import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import mpltern
import pickle
import phase
import os
import sys
import re
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

import argparse
parser = argparse.ArgumentParser(description='This will make images of binodals.')
parser.add_argument('--crits',    metavar='crit',    dest='crits',   type=str, action='store', help='enter address of critical points.'  )
parser.add_argument('--binodals', metavar='bpkl',    dest='bin',     type=str, action='store', help="Name of binodal pkl.", default=None )
parser.add_argument('--mesh',     metavar='mesh',    dest='mesh',     type=str, action='store', help="Name of mesh pkl."                  )
parser.add_argument('--image',    metavar='img',     dest='img',      type=str, action='store', help="Address of image file.",  default=None )
parser.add_argument('--no-rtw',        dest='nrtw',  action='store_true',  help="Dont print out the runtime warning.",           default=False)
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
#########################################

# function from here on out are ravioli-fied versions of what I did in universal.py
def get_crits(P, critpkl):

	if critpkl is None:	
		P.spinodal.obtain_crits()
	else:
		try:
			f = open(critpkl, 'rb')
			crits = pickle.load(f)
			P.spinodal.crits = crits
			f.close()
		except:
			P.spinodal.obtain_crits()
			f = open(critpkl, 'wb')
			pickle.dump(P.spinodal.crits, f)
			f.close()

	P.crits = P.spinodal.crits
	return 
#########################################


if __name__=="__main__":

	start = time.time()

	#####################################################
	crits      = args.crits
	binodals   = args.bin
	images     = args.img
	tern_b       = True
	edges_b      = False
	crits_b      = True

	# get parameters from the crits folder 
	chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", crits)).group(1))
	chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", crits)).group(1))
	chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", crits)).group(1))
	vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", crits)).group(1))
	vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", crits)).group(1))
	vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", crits)).group(1))
	
	params = [chips, chipc, chisc, vs, vc, vp]


	fig = plt.figure(num=0, figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	
	# set up the phase object inputs
	inputs = dict()
	inputs["chi_sc"] = chisc
	inputs["chi_ps"] = chips
	inputs["chi_pc"] = chipc
	inputs["vs"]     = vs
	inputs["vp"]     = vp
	inputs["vc"]     = vc

	print(f"inputs = {inputs}", flush=True)

	# get the phase object
	print(f"Setting up the phase object.")
	P = phase.Phase(inputs)

	get_crits (P, crits)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	f = open(args.bin, 'rb')
	BINODALS = pickle.load(f)
	f.close() 
	
	for idx, H in enumerate(BINODALS["hull_info"]["binodal"]):
		if H[0][0].shape[0] == 0 and H[0][1].shape[0] == 0:
			continue
		elif H[-1] == "two_phase":
			for i in range(0, len(H[0][0]), 100):
				ax.plot([H[0][0][i,0], H[0][1][i,0]], [1-H[0][0][i,0]-H[0][0][i,1], 1-H[0][1][i,0]-H[0][1][i,1]], [H[0][0][i,1], H[0][1][i,1]], c='gold', lw=2, ls='--')
			ax.scatter(H[0][0][:,0], 1-H[0][0][:,0]-H[0][0][:,1], H[0][0][:,1], s=8, c='orangered')
			ax.scatter(H[0][1][:,0], 1-H[0][1][:,0]-H[0][1][:,1], H[0][1][:,1], s=8, c='darkred')

		elif H[-1] == "three_phase":
			ax.plot(np.hstack([H[0][:,0],H[0][0,0]]),\
			np.hstack([1-H[0][:,0]-H[0][:,1], 1-H[0][0,0]-H[0][0,1]]), np.hstack([H[0][:,1], H[0][0,1]]), c='slategray', lw=1)

	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
		
	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

