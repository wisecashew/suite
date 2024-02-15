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
parser.add_argument('--addr-crits',    metavar='acrits',    dest='acrits',   type=str, action='store', help='enter address of critical points.'  )
parser.add_argument('--addr-binodals', metavar='bpkl',      dest='abin',     type=str, action='store', help="Name of binodal pkl.", default=None )
parser.add_argument('--mesh',          metavar='mesh',      dest='mesh',     type=str, action='store', help="Name of mesh pkl."                  )
parser.add_argument('--addr-image',    metavar='img',       dest='img',      type=str, action='store', help="Address of image file.",  default=None )
parser.add_argument('--no-rtw',        dest='nrtw',    action='store_true',  help="Dont print out the runtime warning.",           default=False)
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
	l_crits      = os.listdir(args.acrits)
	l_binodals   = os.listdir(args.abin)
	l_images     = os.listdir(args.img)
	tern_b       = True
	edges_b      = False 
	crits_b      = True

	# get parameters from the crits folder 
	cparams = []
	for num, file in enumerate(l_crits):
		# get the parameters
		chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1))
		chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1))
		chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1))
		vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1))
		cparams.append([chips, chipc, chisc, vs, vc, vp])

	imageparams = []
	for num, file in enumerate(l_images):
		# get the parameters
		chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1))
		chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1))
		chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1))
		vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1))
		imageparams.append([chips, chipc, chisc, vs, vc, vp])

	print(f"Number of binodals is {len(l_binodals)}.", flush=True)

	for num, file in enumerate(l_binodals):
		print(f"At binodal {num}/{len(l_binodals)}...", flush=True)
		# set up the figure 

		# get the parameters
		chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1))
		chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1))
		chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1))
		vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1))
		vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1))

		bparams = [chips, chipc, chisc, vs, vc, vp]
		if bparams in imageparams:
			print(f"Found {bparams} already made.")
			continue

		if bparams in cparams:
			critpkl = l_crits[cparams.index(bparams)]
		else:
			print(f"{bparams} was not found!", flush=True)
			critpkl = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.crits.pkl"

		fig = plt.figure(num=num, figsize=(8,8))
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
		print(f"{args.acrits[:-1]}{critpkl}")
		get_crits (P, f"{args.acrits[:-1]}{critpkl}")

		print(f"Plotting the ternary diagram...", flush=True,end=' ')
		P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
		print(f"done!", flush=True)

		f = open(f"{args.abin[:-1]}{file}", 'rb')
		BINODALS = pickle.load(f)
		f.close() 
		
		for idx, H in enumerate(BINODALS["hull_info"]["binodal"]):
			if H[0][0].shape[0] == 0 and H[0][1].shape[0] == 0:
				continue
			elif H[-1] == "two_phase":
				ax.scatter(H[0][0][:,0], 1-H[0][0][:,0]-H[0][0][:,1], H[0][0][:,1], s=0.5, c='black')
				ax.scatter(H[0][1][:,0], 1-H[0][1][:,0]-H[0][1][:,1], H[0][1][:,1], s=0.5, c='white')
			elif H[-1] == "three_phase":
				ax.plot(np.hstack([H[0][:,0],H[0][0,0]]),\
				np.hstack([1-H[0][:,0]-H[0][:,1], 1-H[0][0,0]-H[0][0,1]]), np.hstack([H[0][:,1], H[0][0,1]]), c='slategray', lw=1)

		fig.savefig(f"{args.img[:-1]}/binimage-chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.png", dpi=1200, bbox_inches="tight")
		
	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

