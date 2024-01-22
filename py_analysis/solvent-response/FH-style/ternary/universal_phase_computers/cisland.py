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
import ternary
import time
import pickle 
import copy
import argparse 
import os 

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
parser.add_argument('--island-stable-pkl',          dest='spkl',       action='store', type=str,  help='File to store stable island information in.')
parser.add_argument('--island-unstable-pkl',        dest='upkl',       action='store', type=str,  help='File to store unstable island information in.')
parser.add_argument('--mesh-pkl',                   dest='mpkl',       action='store', type=str,  help='File to store the mesh in.', default=None)
parser.add_argument('--img-name',                   dest='img',       action='store', type=str,  default="None", help='name of the image to be created (default: all of the inputs in the imagename).')
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
	# phi_s = np.linspace(0.001, 0.999, args.sd)
	# phi_p = np.linspace(0.001, 0.999, args.sd)
	phi_s = np.hstack((np.logspace(-15, -3, 100), np.linspace(1e-3, 0.998, args.sd-200), 1-np.logspace(-3, -15, 100))) # np.logspace(np.log10(0.999), np.log10(1-(1e-15)), 100)))
	phi_p = np.hstack((np.logspace(-15, -3, 100), np.linspace(1e-3, 0.998, args.sd-200), 1-np.logspace(-3, -15, 100))) # np.logspace(np.log10(0.999), np.log10(1-(1e-15)), 100)))

	phi_s, phi_p = np.meshgrid(phi_s, phi_p)

	if args.mpkl is None:
		f = open(f"chips_{chi_ps}-chipc_{chi_pc}-chisc_{chi_sc}-vs_{vs}-vc_{vc}-vp_{vp}.mesh.pkl", 'wb')
		pickle.dump([phi_s, phi_p], f)
		f.close()
	else:
		f = open(args.mpkl, 'wb')
		pickle.dump([phi_s, phi_p], f)
		f.close()

	vals_  = ternary.stab_crit (phi_s, phi_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)
	print(f"vals_.shape = {vals_.shape}", flush=True)
	
	vals_[vals_ >=0]                = 1
	vals_[vals_ < 0]                = 0
	vals_[np.isnan(vals_)]          = 0
	vals_[np.isinf(vals_)]          = 0
	vals_[phi_s+phi_p >= 1]         = 0

	start = time.time()
	

	islands = ternary.find_islands(vals_)
	f = open(args.spkl, 'wb')
	pickle.dump(islands, f)
	f.close()

	# print(f"islands = {islands}")
	stop = time.time()
	print(f"Time of computation = {stop-start} seconds.")
	print(f"# of stable islands = {len(islands)}")

	# repeat the process to get the number of unstable islands
	# phi_s = np.linspace(0.001, 0.999, args.sd)
	# phi_p = np.linspace(0.001, 0.999, args.sd)
	phi_s = np.hstack((np.logspace(-15, -3, 100), np.linspace(1e-3, 0.998, args.sd-200), 1-np.logspace(-3, -15, 100))) # np.logspace(np.log10(0.999), np.log10(1-(1e-15)), 100)))
	phi_p = np.hstack((np.logspace(-15, -3, 100), np.linspace(1e-3, 0.998, args.sd-200), 1-np.logspace(-3, -15, 100))) # np.logspace(np.log10(0.999), np.log10(1-(1e-15)), 100)))
	
	phi_s, phi_p = np.meshgrid(phi_s, phi_p)

	vals_  = ternary.stab_crit (phi_s, phi_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)
	print(f"vals_.shape = {vals_.shape}", flush=True)
	vals_[vals_ >=0]                 = 0
	vals_[np.isnan(vals_)]           = 0
	vals_[np.isinf(vals_)]           = 0
	vals_[phi_s+phi_p >= 1]          = 0
	vals_[vals_ < 0]                 = 1

	start = time.time()

	islands = ternary.find_islands(vals_)
	f = open(args.upkl, 'wb')
	pickle.dump(islands, f)
	f.close()

	if len(islands) == 0:
		os.rename(args.spkl, args.spkl+"_TODEL")
		os.rename(args.upkl, args.upkl+"_TODEL")

	stop = time.time()
	print(f"Time of computation = {stop-start} seconds.")
	print(f"# of unstable islands = {len(islands)}")

	###################################
	# islands have been dumped out
	###################################

	to_keep = ~np.isnan(ternary.stab_crit(p_s, p_p, vs, vc, vp, chi_ps, chi_pc, chi_sc))
	vals = vals [to_keep]
	p_s  = p_s  [to_keep]
	p_p  = p_p  [to_keep]

	if len(vals) == 0:
		print (f"There will be no critical points and no spinodal region.")

	start_color = 'indianred'
	end_color   = 'steelblue'

	# Number of steps in the gradient
	cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [start_color, end_color])
	vmax = +1
	vmin = -1
	norm = colors.Normalize(vmin=vmin, vmax=vmax)
	vals[vals>=0] = 1
	vals[vals<0]  = 0
	cols = cmap (norm (vals))

	if len(vals) == 0:
		print (f"There will be no critical points and no spinodal region.", flush=True)

	if (vals > 0).all():
		print(f"There are no unstable regions.", flush=True)
	elif (vals>=0).any() and (vals<0).any():
		print(f"There are stable and unstable regions")

	# get all the crit points
	if args.crits:
		roots_up, roots_down = ternary.find_crit_point(vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s)

		# put them all together
		crits      = np.vstack ((roots_up, roots_down))

		# get rid ofthe redundant ones
		threshold  = 1e-6
		crits, keep      = ternary.remove_close_rows (crits, threshold)

		print (f"cleaned_crits = \n{crits}", flush=True)

	else:
		crits = None

	# Plot the points
	p_c = 1 - p_s - p_p

	# plot the thing
	# print(f"ax={ax}, args.ternary={args.ternary}, args.edges={args.edges}, args.pc={args.pc}, crits={crits}, chi_ps, chi_pc, chi_sc={chi_ps, chi_pc, chi_sc}, cols={cols[0:5]}, root_up_p, root_lo_p={root_up_p,root_lo_p}, root_up_s, root_lo_s = {root_up_s, root_lo_s}", flush=True)
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



