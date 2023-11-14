import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import tangent
import mpltern
import phase

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--plot-edges', dest='pe',    action='store_true', help='plot the edges of the spinodal.')
parser.add_argument('--plot-crits', dest='pc',    action='store_true', help='plot the critical points.')
parser.add_argument('--no-rtw',     dest='nrtw',  action='store_true',  default=False, help="Dont print out the runtime warning.")
parser.add_argument('--img',        dest='img',   action='store',      type=str, default="None", help="name of image to be created. (default: blastradius).")
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

if __name__=="__main__":

	start = time.time()

	#####################################################
	fig = plt.figure(num=0, figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	########################################################
	# set up the inputs
	# set up the chi values
	inputs = dict()
	inputs["chi_sc"] = args.chi_sc
	inputs["chi_ps"] = args.chi_ps
	inputs["chi_pc"] = args.chi_pc
	inputs["vs"] = args.vs
	inputs["vp"] = args.vp
	inputs["vc"] = args.vc
	tern_b  = True
	edges_b = args.pe
	crits_b = args.pc

	# set up objects and crit points
	print(f"Setting up objects...", flush=True, end=' ')
	P = phase.Phase(inputs)
	print("done!", flush=True)
	P.spinodal.obtain_crits()
	P.crits = P.spinodal.crits
	print(f"crits = {P.crits}")
	print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	# start charting out the curve
	
	for idx in range(len(P.crits)):
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_p1, phi_p2 = P.sym_mu_ps.binodal_run_in_s2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='#C8442F', s=1)
		ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='#2FB3C8', s=1)

		phi_s1, phi_s2, phi_p1, phi_p2 = P.sym_mu_ps.binodal_run_in_p2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='#C8442F', s=1)
		ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='#2FB3C8', s=1)

		print("Plotted out!", flush=True)

	for idx in range(len(P.crits)):
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_c1, phi_c2 = P.sym_mu_sc.binodal_run_in_s2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='#346ECB',   s=1)
		ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='#CB9134', s=1)

		phi_s1, phi_s2, phi_c1, phi_c2 = P.sym_mu_sc.binodal_run_in_c2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='#346ECB',   s=1)
		ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='#CB9134', s=1)

		print("Plotted out!", flush=True)

	for idx in range(len(P.crits)):
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_p1, phi_p2, phi_c1, phi_c2 = P.sym_mu_pc.binodal_run_in_p2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(1-phi_c1-phi_p1, phi_c1, phi_p1, c='#C8442F',   s=1)
		ax.scatter(1-phi_c2-phi_p2, phi_c2, phi_p2, c='#2FB3C8', s=1)

		phi_p1, phi_p2, phi_c1, phi_c2 = P.sym_mu_pc.binodal_run_in_c2(P.crits[idx])

		# start plotting out the curve
		ax.scatter(1-phi_c1-phi_p1, phi_c1, phi_p1, c='#C8442F',   s=1)
		ax.scatter(1-phi_c2-phi_p2, phi_c2, phi_p2, c='#2FB3C8', s=1)

		print("Plotted out!", flush=True)

	if len(P.crits) == 1:
		along_normal = True
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_p1, phi_p2 = P.sym_mu_ps.binodal_run_in_s2(P.crits[idx], along_normal)

		# start plotting out the curve
		ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='#C8442F',   s=1)
		ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='#2FB3C8', s=1)

		phi_s1, phi_s2, phi_p1, phi_p2 = P.sym_mu_ps.binodal_run_in_p2(P.crits[idx], along_normal)

		# start plotting out the curve
		ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='#C8442F', s=1)
		ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='#2FB3C8', s=1)

		print("Plotted out!", flush=True)

	# create the image
	print("Making image...")
	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		fig.savefig (img_name, dpi=1200, bbox_inches="tight")
	else:
		fig.savefig (f"bin_tern-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	print(f"done!", flush=True)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

