import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
from scipy.spatial.distance import cdist
import sys
import argparse
import linecache
import mpltern
np.set_printoptions(threshold=sys.maxsize)
import warnings 
import tangent
import ternary
import time
from sklearn.cluster import DBSCAN
import pickle 
import copy
import argparse 

parser = argparse.ArgumentParser(description="Locate the /c/ritical points on the /spin/odal diagram. This program will create one plot, and you can customize what you want on the plot.")
parser.add_argument('--chisc',  metavar='chi_sc',  dest='chi_sc',  type=float,   action='store', help='enter A-C exchange parameter.' )
parser.add_argument('--chips',  metavar='chi_ps',  dest='chi_ps',  type=float,   action='store', help='enter A-B exchange parameter.' )
parser.add_argument('--chipc',  metavar='chi_pc',  dest='chi_pc',  type=float,   action='store', help='enter B-C exchange parameter.' )
parser.add_argument('-vs',      metavar='vs',      dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',      dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',      dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--search-density',      metavar='sd',      dest='sd',      type=int,   action='store', help='density at which the search should be performed for stable regions.',default=100)
parser.add_argument('--dont-calc-crits',     dest='crits',     action='store_false', default=True,  help='Put this in to make sure critical points are not calculated.')
parser.add_argument('--ternary',             dest='ternary',   action='store_true',  default=False, help='make the output a ternary plot.')
parser.add_argument('--plot-crits',          dest='pc',        action='store_true',  default=False, help='plot critical points.')
parser.add_argument('--draw-edges-spinodal', dest='edges',     action='store_true',  default=False, help='draw the edges of the spinodal.')
parser.add_argument('--no-rtw',              dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--tang_norm',           dest='tang_norm', action='store_true',  default=False, help='draw normal and tangent at critical point.')
parser.add_argument('--pkl',                 dest='pkl',       action='store', type=str,  help='draw normal and tangent at critical point.')
parser.add_argument('--img-name',            dest='img',       action='store', type=str,  default="None", help='name of the image to be created (default: all of the inputs in the imagename).')
args = parser.parse_args()

def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f"beep.\n"
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format


if __name__=="__main__":

	print ("Revving up the program...", flush=True)

	####################################################
	# set up the inputs
	lsize = 3
	font = {'color':  'black',
		'weight': 'normal',
		'size': lsize}

	fig = plt.figure(figsize=(8,8))
	if args.ternary:
		ax = fig.add_subplot(projection="ternary")
	else:
		ax = plt.axes()

	chi_sc = args.chi_sc
	chi_ps = args.chi_ps
	chi_pc = args.chi_pc
	vs     = args.vs
	vc     = args.vc
	vp     = args.vp
	####################################################


	# solution in terms of phi_s
	discriminant_s = lambda phi_s: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*\
	vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc))+(vp - 2*phi_s*vp *vs *chi_ps + \
	vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 \
	- 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))))**2

	denom_s    = lambda phi_s: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))
	prefac_s   = lambda phi_s: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp *\
	(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc)))

	root_up_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) + np.sqrt(discriminant_s(phi_s)))
	root_lo_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) - np.sqrt(discriminant_s(phi_s)))


    # solution in terms of phi_p
	discriminant_p = lambda phi_p: -4*vc*vs*(phi_p*vp+(-1+phi_p)*vc*(-1+2*phi_p*vp*chi_pc))*(2*chi_sc+phi_p*vp*\
	(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc)))+(vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*\
	vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2\
	-2*chi_pc*(1+vs*(chi_ps+chi_sc))))))**2

	denom_p        = lambda phi_p: 1/(-2*vc*vs*(2*chi_sc+phi_p*vp*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))))
	prefac_p       = lambda phi_p: vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-\
	2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2-2*chi_pc*(1+vs*(chi_ps+chi_sc)))))

	root_up_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)+np.sqrt(discriminant_p(phi_p)))
	root_lo_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)-np.sqrt(discriminant_p(phi_p)))


	p_s_space = np.arange (0.001, 1-0.001, 0.001)
	p_s = np.repeat (p_s_space, len(p_s_space))

	p_p = np.zeros (p_s.shape)
	for i in range (len(p_s_space)):
		p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

	vals = ternary.stab_crit(p_s, p_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)


	###############################
	# stuff to make the islands
	###############################
	phi_s = np.linspace(0.001, 0.999, args.sd)
	phi_p = np.linspace(0.001, 0.999, args.sd)
	
	phi_s, phi_p = np.meshgrid(phi_s, phi_p)

	vals_  = ternary.stab_crit (phi_s, phi_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)
	print(f"vals_.shape = {vals_.shape}", flush=True)
	# print(f"vals_ = {vals_}")
	vals_[vals_ >=0]            = 1
	vals_[vals_ < 0]            = 0
	vals_[np.isnan(vals_)]      = 0
	vals_[np.isinf(vals_)]      = 0
	vals_[phi_s+phi_p >= 0.9999999] = 0

	start = time.time()

	# print(f"vals = {vals}")
	islands = ternary.find_islands(vals_)

	print(f"islands = {islands}")
	stop = time.time()
	print(f"Time of computation = {stop-start} seconds.")
	print(f"# of islands = {len(islands)}")
	f = open(args.pkl, 'wb')
	pickle.dump(islands, f)
	f.close()
	###################################
	# islands have been dumped out
	###################################

	to_keep = ~np.isnan(ternary.stab_crit(p_s, p_p, vs, vc, vp, chi_ps, chi_pc, chi_sc))
	vals = vals [to_keep]
	p_s  = p_s  [to_keep]
	p_p  = p_p  [to_keep]

	if len(vals) == 0:
		print (f"There will be no critical points and no spinodal region.")

	vmax = np.max (vals)
	vmin = np.min (vals)
	norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
	cols = cm.bwr (norm (vals))

	if np.sign (vmax) == np.sign (vmin):
		if np.sign (vmax) >=0:
			vmin = -vmax
			print (f"There is no unstable region.", flush=True)
		else:
			vmax = -vmin
			print ("There is mostly unstable region.", flush=True)

	else:
		print ("there exist unstable regions.", flush=True)

	# get all the crit points
	if args.crits:
		roots_up, roots_down = ternary.find_crit_point(vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s)

		# put them all together
		crits      = np.vstack ((roots_up, roots_down))

		# get rid ofthe redundant ones
		threshold  = 1e-6
		crits      = ternary.remove_close_rows (crits, threshold)

		print (f"cleaned_crits = \n{crits}", flush=True)

	else:
		crits = None

	# Plot the points
	p_c = 1 - p_s - p_p

	# plot the thing
	spinodal_phi_s, spinodal_phi_p = ternary.plot(ax, args.ternary, args.edges, args.pc, crits, chi_ps, chi_pc, chi_sc, p_s, p_p, cols, root_up_p, root_lo_p, root_up_s, root_lo_s)
	pts = np.array([spinodal_phi_s, spinodal_phi_p]).T

	# plot the tangent and normal
	ternary.add_tang_norm(ax, args.tang_norm, args.ternary, args.pc, crits, vs, vc, vp, chi_pc, chi_ps, chi_sc, root_up_s, root_lo_s)

	# add the plot embellishments
	ternary.embelish(ax, args.ternary)

	# plot the grid on the axis
	ax.grid()

	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		plt.savefig (img_name, dpi=1200, bbox_inches="tight")
	elif args.ternary:
		plt.savefig (f"signs_tern-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)
	else:
		plt.savefig (f"signs_reg-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)

	print ("Completed heat map computation.", flush=True)



