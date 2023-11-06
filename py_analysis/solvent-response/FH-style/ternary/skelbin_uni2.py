import numpy as np
import numba as nb
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
from scipy.optimize import root
import scipy.spatial.distance
from scipy.spatial.distance import cdist
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import sys
import os
import argparse
import time
import warnings
import linecache
import itertools
import ternary
import tangent

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--no-rtw', dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--mesh',   metavar='mesh',   dest='mesh', type=int, action='store', help='enter mesh fineness.')
parser.add_argument('--skelfile', dest='skelfile', type=str, action='store', help="name of file where we are going to dump the skeleton of the binodal (the default contains all the information provided above).", default="None")
parser.add_argument('--out-spinodal', dest='ons', action='store_true', default=False, help="Keep fractions that are outside the spinodal (default: False).")
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    if args.nrtw:
        return f"beep.\n"
    else:
        return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

#########################################
def perform_sweep (phi_b, mesh, chi_ps, chi_pc, chi_sc, crit_point, phi_a_edge, phi_b_edge, center):

	print (f"pid = {os.getpid()}.", flush=True)
	phi_b = np.repeat (phi_b, mesh)
	phi_a = np.zeros  (phi_b.shape)

	for i in range (mesh):
		upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
		phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

	# only keep stuff which is outside the spinodal
	if args.ons:
		to_keep = ternary.stab_crit (phi_a, phi_b, chi_ps, chi_pc, chi_sc) >= 0
	else:
		to_keep = np.array([True]*len(phi_a))

# define the fractions along which 
	phi_b   = phi_b [to_keep]
	phi_a   = phi_a [to_keep]

	phi_b   = np.hstack((phi_b, phi_b_edge))
	phi_a   = np.hstack((phi_a, phi_a_edge))

# now start splitting up phi_a, phi_b
	phis    = np.vstack((phi_a, phi_b)).T

# define the axis along which to make the split
	central_axis   = (crit_point-center) / np.linalg.norm (crit_point-center)

	# find those ABOVE axis
	direction      = (phis - center) / np.linalg.norm(phis-center, axis=1)[:, np.newaxis]
	clock          = np.cross (central_axis, direction)
	phi_upper      = phis[clock > 0]
	phi_lower      = phis[clock < 0]

	# now split phi_a, phi_b such that they are on either side of the critical line
	chem_pot_a_upper = mu_a (phi_upper[:,0], phi_upper[:,1])
	chem_pot_b_upper = mu_b (phi_upper[:,0], phi_upper[:,1])
	chem_pot_c_upper = mu_c (phi_upper[:,0], phi_upper[:,1])

	masks = np.isinf(chem_pot_a_upper) | np.isnan(chem_pot_a_upper) | np.isinf(chem_pot_b_upper) |\
	 np.isnan(chem_pot_b_upper) | np.isinf (chem_pot_c_upper) | np.isnan (chem_pot_c_upper)

	chem_pot_a_upper = chem_pot_a_upper [~masks]
	chem_pot_b_upper = chem_pot_b_upper [~masks]
	chem_pot_c_upper = chem_pot_c_upper [~masks]

	# get only the relevant phi_uppers
	phi_upper = phi_upper[~masks]

	# calculate the chemical potentials of the lower half
	chem_pot_a_lower = mu_a (phi_lower[:,0], phi_lower[:,1])
	chem_pot_b_lower = mu_b (phi_lower[:,0], phi_lower[:,1])
	chem_pot_c_lower = mu_c (phi_lower[:,0], phi_lower[:,1])

	# create a mask for all the good lower chemical potentials
	masks = np.isinf(chem_pot_a_lower) | np.isnan(chem_pot_a_lower) | np.isinf(chem_pot_b_lower) |\
	np.isnan(chem_pot_b_lower) | np.isinf (chem_pot_c_lower) | np.isnan (chem_pot_c_lower)

	# get only the relevant phi_lowers
	chem_pot_a_lower = chem_pot_a_lower [~masks]
	chem_pot_b_lower = chem_pot_b_lower [~masks]
	chem_pot_c_lower = chem_pot_c_lower [~masks]
	phi_lower = phi_lower[~masks]

	# get an array of these values
	mu_upper   = np.array ([chem_pot_a_upper, chem_pot_b_upper, chem_pot_c_upper]).T
	mu_lower   = np.array ([chem_pot_a_lower, chem_pot_b_lower, chem_pot_c_lower]).T

	# calculate all the distances
	distances       = np.linalg.norm (mu_upper[:, np.newaxis] - mu_lower, axis=2)
	print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
	closest_indices = np.argmin(distances, axis=1)

	# partition, partition... 
	phi_lower       = phi_lower[closest_indices]
	mu_lower        = mu_lower [closest_indices]
	min_distances   = distances[np.arange(len(mu_upper)), closest_indices]

	return [phi_upper, phi_lower, min_distances]

#########################################

if __name__=="__main__":

	start = time.time()

	############################################################################################################################

	chi_sc = args.chi_sc
	chi_ps = args.chi_ps
	chi_pc = args.chi_pc

	vs = args.vs
	vp = args.vp
	vc = args.vc

	############################################################################################################################
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

	############################################################################################################################
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

	############################################################################################################################
	# find the edges of the spinodal curve
	meshsize            = 1000
	phi_s               = np.linspace (0.001, 1-0.001, meshsize*10)

	r1 = root_up_s (phi_s)
	r2 = root_lo_s (phi_s)

	to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
	r1 = r1[to_keep_1]

	to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
	r2 = r2[to_keep_2]

	phi_a_edge_1 = phi_s[to_keep_1]
	phi_a_edge_2 = phi_s[to_keep_2]

	phi_a_edge = np.hstack((phi_a_edge_1,phi_a_edge_2))
	phi_b_edge = np.hstack((r1,r2))
	if len(phi_a_edge)==0:
		print("There is no critical region. No binodals will be found. Exiting...")
		exit()

	############################################################################################################################
	# get the critical points
	roots_up, roots_down = ternary.find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s)
	crits = np.vstack((roots_up, roots_down))
	crits = ternary.remove_close_rows (crits, 1e-6)

	# FIND PHI_B GIVEN PHI_A
	mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - vs/vp * phi_b - vs/vc * (1-phi_a-phi_b) + vs * (phi_b**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_sc + phi_b * (1-phi_a-phi_b) * (chi_ps + chi_sc - chi_pc) ) 
	mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vp/vs * phi_a - vp/vc * (1-phi_a-phi_b) + vp * (phi_a**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_pc + phi_a * (1-phi_a-phi_b) * (chi_ps + chi_pc - chi_sc) )
	mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/vs * phi_a - vc/vp * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )

	mesh  = args.mesh

	if args.skelfile == "None":
		skelfile = f"bin-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.skelfile"
	else:
		skelfile = args.skelfile
	f = open (skelfile, 'w')
	f.write ("{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}\n"\
	.format ("dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2")) 

	phi_b_spin_min = 0.01
	phi_b_spin_max = 0.8

	phi_b_list = [np.linspace(phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1, mesh), \
		np.linspace(phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3, mesh), \
		np.linspace(phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5, mesh), \
		np.linspace(phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9, mesh), \
		np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh), \
		np.logspace(-20, np.log10(0.9), mesh)]

	print (f"crits = \n{crits}", flush=True)
	for crit in crits:
		print (f"@ crit point = {crit}...", flush=True)
		tangent_to_crit = tangent.tangent2(vs, vc, vp, crit[0], crit[1], chi_ps, chi_pc, chi_sc, root_up_s, root_lo_s)
		normal_slope    = -1/tangent_to_crit
		if len(crits) == 1:
			center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit
		elif len(crits) == 2:
			center  = np.mean(crits, axis=0)
		else:
			center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit

		f.write  ("{:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f}\n"\
		.format(0, crit[0], crit[1], 1-crit[0]-crit[1], crit[0], crit[1], 1-crit[0]-crit[1]))


		print (f"# of iterations = {len(phi_b_list)}...", flush=True)
		for idx,b_list in enumerate(phi_b_list):
			print (f"\tAt iteration {idx}...", flush=True)
			results = perform_sweep (b_list, mesh, chi_ps, chi_pc, chi_sc, crit, phi_a_edge, phi_b_edge, center)

			if len(results[0]) == 0:
				print (f"\tNo critical condition in process {idx}.", flush=True)
				continue
			for i in range(len (results[0]) ):
				f.write  ("{:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f}\n"\
		.format(results[2][i], results[0][i,0], results[0][i,1], 1-results[0][i,0]-results[0][i,1], results[1][i,0], results[1][i,1], 1-results[1][i,0]-results[1][i,1]))

	f.close ()

	stop = time.time ()
	print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)


