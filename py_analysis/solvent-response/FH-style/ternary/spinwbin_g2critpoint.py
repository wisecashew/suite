import numpy as np
import pandas as pd
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
from scipy.spatial import ConvexHull
import sys
import argparse
import time
import multiprocessing as mp
import mpltern
import os
import itertools
import warnings
import linecache
import tangent
import ternary

np.set_printoptions(threshold=sys.maxsize)
# os.system("taskset -p 0xfffff %d" % os.getpid())
# os.environ['MKL_NUM_THREADS'] = '1'
# os.environ['NUMEXPR_NUM_THREADS'] = '1'
# os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush()


import argparse
parser = argparse.ArgumentParser(description="Create a ternary /spin/odal and /bin/odal diagram.")
parser.add_argument('--chisc',           metavar='chi_sc',    dest='chi_sc',       type=float,     action='store',             help='enter S-C exchange parameter.')
parser.add_argument('--chips',           metavar='chi_ps',    dest='chi_ps',       type=float,     action='store',             help='enter P-S exchange parameter.')
parser.add_argument('--chipc',           metavar='chi_pc',    dest='chi_pc',       type=float,     action='store',             help='enter P-C exchange parameter.')
parser.add_argument('-vs',               metavar='vs',        dest='vs',           type=float,     action='store',             help='specific volume of solvent.')
parser.add_argument('-vc',               metavar='vc',        dest='vc',           type=float,     action='store',             help='specific volume of cosolvent.')
parser.add_argument('-vp',               metavar='vp',        dest='vp',           type=float,     action='store',             help='specific volume of polymer.')
parser.add_argument('--nadded-rows',     metavar='nar',       dest='nar',          type=int,       action='store',             help='number of rows to add while refining.', default=10)
parser.add_argument('--refine-count',    metavar='rc',        dest='rc',           type=int,       action='store',             help='number of times the refiner should go off (default: 10).', default=10)
parser.add_argument('--refine-by-size',  metavar='rbs',       dest='rbs',          type=int,       action='store',             help='Refine the binodal by size. (default: -1)', default=-1)
parser.add_argument('--diff-count',      metavar='dc',        dest='dc',           type=int,       action='store',             help='number of times the difference filler should go off (default: 10).', default=10)
parser.add_argument('--tieline-density', metavar='td',        dest='td',           type=int,       action='store',             help='Plug in a density at which you want tielines (default: 50).', default=50)
parser.add_argument('--dumpfile',        dest='dumpfile',     type=str,            action='store', help="name of file where the skeleton was dumped.")
parser.add_argument('--bin-boundary',    dest='boundary',     type=str,            action='store', help="name of file where you will dump out your solution for the binodal (default name holds information about all inputs).", default="None")
parser.add_argument('--ternary',         action='store_true', default=False,       help='make the output a ternary plot.')
parser.add_argument('--only-one-crit',   dest='ooc',          action='store_true', default=False,  help='only go through one crit point. Useful only while debugging, probably. Avoid using otherwise! (default: False)')
parser.add_argument('--tielines',        dest='tl',           action='store_true', default=False,  help="Option to include if you want to see all tie-lines.")
parser.add_argument('--out-spinodal',    dest='ous',          action='store_true', default=False,  help="Only consider compositions that are outside of the spinodal (default: False).")
parser.add_argument('--no-rtw',          dest='nrtw',         action='store_true', default=False,  help="Don't print out the runtime warning.")
parser.add_argument('--normal-cleave',   dest='normal',       action='store_true', default=False,  help="Cleave along the normal (default: cleave along the center of mass of critical points).")
parser.add_argument('--image',           dest='img',          type=str,            action='store', help="name of image generated (default name holds information about all inputs).", default="None")
args = parser.parse_args()

if args.boundary == "None":
	boundaryfile = f"vs_{args.vs}-vc_{args.vc}-vp_{args.vp}-chisc_{args.chi_sc}-chips_{args.chi_ps}-chipc_{args.chi_pc}.binodal"
else:
	boundaryfile = args.boundary

###############

######
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f"beep.\n"
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

################################

###############################

def add_interpolated_rows (original_array, M):
	result_array = []

	for i in range(len(original_array) - 1):
		start_row = original_array[i]
		end_row = original_array[i + 1]
		step = (end_row - start_row) / (M + 1)

		interpolated_rows = [start_row + step * j for j in range(1, M + 1)]
		result_array.extend([start_row] + interpolated_rows)

	result_array.append(original_array[-1])
	return np.array(result_array)

##############################

def add_rows_between_largest_gap (array, M):

	diff        = np.diff (array, axis=0)
	distances   = np.linalg.norm (diff, axis=1)
	max_dists   = np.argmax (distances)
	insert_rows = np.linspace (array[max_dists], array[max_dists+1], M)
	array       = np.insert (array, max_dists+1, insert_rows[1:-1], axis=0)

	return array, max_dists

###############################

def add_rows_at_index (array, idx, M):

	insert_rows = np.linspace (array[idx], array[idx+1], M)
	array       = np.insert (array, idx+1, insert_rows[1:-1],axis=0)

	return array

###############################

def max_dists_on_binodal (top_half, bottom_half):

	diff_top         = np.diff        (top_half, axis=0)
	distances        = np.linalg.norm (diff_top, axis=1)
	max_dists_top    = np.max         (distances)

	diff_bottom      = np.diff        (bottom_half, axis=0)
	distances        = np.linalg.norm (diff_bottom, axis=1)
	max_dists_bottom = np.max         (distances)

	max_dist = np.max ([max_dists_top, max_dists_bottom])

	return max_dist

######################################

def refined_binodal_v4 (side_1, side_2, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc):

	side_1x, m1 = add_rows_between_largest_gap (side_1, nadded_rows)
	side_2x     = add_rows_at_index (side_2, m1, nadded_rows)

	print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v4...", flush=True)
# print (f"m1 = {m1}, m2 = {m2}.")
	print (f"from {side_1x[m1+1]} to {side_1x[m1+nadded_rows-1]}", flush=True)

	add_counter = 0
	for idx, pt in enumerate (side_1x[m1+1:m1+nadded_rows-1]):
		if idx%25==0: print (f"idx = {idx} @ x,y = {pt[0],pt[1]}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root_store = []
		dist_store = []

		for tidx, tpt in enumerate (side_2x[m1+1+idx-nadded_rows:m1+1+idx+nadded_rows]):
			root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]])
			# if the roots are "bad" roots, just write them out as bad
			if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
				continue

			else:
				fa = [pt[0], root[1]]
				fb = [root[0], root[2]]
				fc = [1-pt[0]-root[0], 1-root[1]-root[2]]
				p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						root_store.append((p1,p2))
					else:
						root_store.append((p2,p1))
					dist_store.append (np.linalg.norm (p1-p2))

		# choose the root that is furthest away
		try:
			best_root  = np.argmax (dist_store)
			root_combo = root_store[best_root]
			side_1 = np.insert (side_1, m1+1+add_counter, root_combo[0], axis=0)
			side_2 = np.insert (side_2, m1+1+add_counter, root_combo[1], axis=0)
			add_counter += 1
			# side_1[m1+1+idx] = root_combo[0]
			# side_2[m1+1+idx] = root_combo[1]

		except:
			print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

	return [side_1, side_2]

#####################################

def refined_binodal_v6 (side_1, side_2, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc):

	side_2x, m2 = add_rows_between_largest_gap (side_2, nadded_rows)
	side_1x     = add_rows_at_index (side_1, m2, nadded_rows)

	print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v6...", flush=True)
# print (f"m1 = {m1}, m2 = {m2}.")
	print (f"from {side_2x[m2+1]} to {side_2x[m2+nadded_rows-1]}", flush=True)

	add_counter = 0
	for idx, pt in enumerate (side_2x[m2+1:m2+nadded_rows-1]):
		if idx%25==0: print (f"idx = {idx} @ x,y = {pt[0],pt[1]}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root_store = []
		dist_store = []

		for tidx, tpt in enumerate (side_1x[m2+1+idx-nadded_rows:m2+1+idx+nadded_rows]):
			root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]])

			# if the roots are "bad" roots, just write them out as bad
			if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
				continue

			else:
				fa = [pt[0], root[1]]
				fb = [root[0], root[2]]
				fc = [1-pt[0]-root[0], 1-root[1]-root[2]]
				p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						root_store.append((p1,p2))
					else:
						root_store.append((p2,p1))
					dist_store.append (np.linalg.norm (p1-p2))

		# choose the root that is furthest away
		try:
			best_root  = np.argmax (dist_store)
			root_combo = root_store[best_root]
			side_1 = np.insert (side_1, m1+1+add_counter, root_combo[1], axis=0)
			side_2 = np.insert (side_2, m1+1+add_counter, root_combo[0], axis=0)
			add_counter += 1

		except:
			print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

	return [side_1, side_2]
######################################

def refined_binodal_v7 (side_1, side_2, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc):

	side_2x, m2 = add_rows_between_largest_gap (side_2, nadded_rows)
	side_1x     = add_rows_at_index (side_1, m2, nadded_rows)

	print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v7...", flush=True)
# print (f"m1 = {m1}.")
	print (f"from {side_2x[m2+1]} to {side_2x[m2+nadded_rows-1]}", flush=True)
	add_counter = 0

	for idx, pt in enumerate (side_2x[m2+1:m2+nadded_rows-1]):
		if idx%25==0: print (f"idx = {idx} @ x, y = {pt[0],pt[1]}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(phi[0], pt[1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], pt[1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], pt[1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root_store = []
		dist_store = []

		for tidx, tpt in enumerate (side_1x[m2+1+idx-2*nadded_rows:m2+1+idx+2*nadded_rows]):
			root = fsolve (mu_equations, [pt[0], tpt[0], tpt[1]])

			# if the roots are "bad" roots, just write them out as bad
			if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
				continue

			else:
				fa = [root[0], root[1]]
				fb = [pt[1]  , root[2]]
				fc = [1-root[0]-pt[1], 1-root[1]-root[2]]
				p1 = np.array([root[0], pt[1], 1-root[0]-pt[1]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue
				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue
				else: # if the roots are basically the same point, write them out as bad 
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						root_store.append((p1,p2))
					else:
						root_store.append((p2,p1))
					dist_store.append (np.linalg.norm (p1-p2))

		# choose the root that is furthest away
		try:
			best_root = np.argmax (dist_store)
			root_combo = root_store[best_root]
			side_1 = np.insert (side_1, m2+1+add_counter, root_combo[1], axis=0)
			side_2 = np.insert (side_2, m2+1+add_counter, root_combo[0], axis=0)
			add_counter += 1
			# side_1[m1+1+idx] = root_combo[0]
			# side_2[m1+1+idx] = root_combo[1]

		except:
			print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

	return [side_1, side_2]

######################################


######################################

def refined_binodal_v5 (side_1, side_2, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc):

	side_1x, m1 = add_rows_between_largest_gap (side_1, nadded_rows)
	side_2x     = add_rows_at_index (side_2, m1, nadded_rows)

	print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v5...", flush=True)
# print (f"m1 = {m1}.")
	print (f"from {side_1x[m1+1]} to {side_1x[m1+nadded_rows-1]}", flush=True)
	add_counter = 0

	for idx, pt in enumerate (side_1x[m1+1:m1+nadded_rows-1]):
		if idx%25==0: print (f"idx = {idx} @ x, y = {pt[0],pt[1]}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(phi[0], pt[1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], pt[1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], pt[1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root_store = []
		dist_store = []

		for tidx, tpt in enumerate (side_2x[m1+1+idx-2*nadded_rows:m1+1+idx+2*nadded_rows]):
			root = fsolve (mu_equations, [pt[0], tpt[0], tpt[1]])

			# if the roots are "bad" roots, just write them out as bad
			if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
				continue

			else:
				fa = [root[0], root[1]]
				fb = [pt[1]  , root[2]]
				fc = [1-root[0]-pt[1], 1-root[1]-root[2]]
				p1 = np.array([root[0], pt[1], 1-root[0]-pt[1]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else: # if the roots are basically the same point, write them out as bad 
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						root_store.append((p1,p2))
					else:
						root_store.append((p2,p1))
					dist_store.append (np.linalg.norm (p1-p2))

		# choose the root that is furthest away
		try:
			best_root = np.argmax (dist_store)
			root_combo = root_store[best_root]
			side_1 = np.insert (side_1, m1+1+add_counter, root_combo[0], axis=0)
			side_2 = np.insert (side_2, m1+1+add_counter, root_combo[1], axis=0)
			add_counter += 1
			# side_1[m1+1+idx] = root_combo[0]
			# side_2[m1+1+idx] = root_combo[1]

		except:
			print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

	return [side_1, side_2]

######################################

def root_finder_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string, iterr, chi_ps, chi_pc, chi_sc):

	if scale_string == "phi_a":
		scale_a = 1e+6*(10**iterr)
		scale_b = 1
	elif scale_string == "phi_b":
		scale_a = 1
		scale_b = 1e+6*(10**iterr)
	else:
		print (f"Bad string provided: {scale_string}.")
		exit()

	direction     = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_upper   = np.arccos(np.dot(direction, central_axis))
	theta1        = theta_upper [max_ind]
	theta2        = theta_upper [max_ind+1]

	unsolved_upper_bin = binodal_upper[bad_idx][:, 0:2]
	direction_binodal  = (unsolved_upper_bin[:,0:2] - center) / np.linalg.norm(unsolved_upper_bin[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_binodal      = np.arccos(np.dot(direction_binodal, central_axis))
	to_keep            = (theta_binodal > theta1) & (theta_binodal < theta2)

	unsolved_upper_bin = unsolved_upper_bin [to_keep]

	sol_bin_up    = np.empty((0,3))
	sol_bin_down  = np.empty((0,3))

	lower_guesses = np.linspace (sol_lower[max_ind], sol_lower[max_ind+1], 100)

	for idx in range (0, len(unsolved_upper_bin), 100):
		print (f"@ idx = {idx}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(phi[0], unsolved_upper_bin[idx][1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], unsolved_upper_bin[idx][1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], unsolved_upper_bin[idx][1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for iidx in range(len(lower_guesses)):
			root = fsolve (mu_equations, [unsolved_upper_bin[idx][0], lower_guesses[iidx][0]/scale_a, lower_guesses[iidx][1]/scale_b])

			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue

			else:
				p1 = np.array([root[0], unsolved_upper_bin[idx][1], 1-root[0]-unsolved_upper_bin[idx][1]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.linalg.norm(p1-p2) < 1e-6:
					continue

				elif np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)

					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up,  p1))
						sol_bin_down = np.vstack((sol_bin_down,p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up,  p2))
						sol_bin_down = np.vstack((sol_bin_down,p1))
					break

	for idx in range (0, len(unsolved_upper_bin), 100):
		print (f"@ idx = {idx}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(unsolved_upper_bin[idx][0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(unsolved_upper_bin[idx][0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(unsolved_upper_bin[idx][0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for iidx in range(len(lower_guesses)):
			root = fsolve (mu_equations, [unsolved_upper_bin[idx][1], lower_guesses[iidx][0]/scale_a, lower_guesses[iidx][1]/scale_b])

			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue

			else:
				p1 = np.array([unsolved_upper_bin[idx][0], root[0], 1-root[0]-unsolved_upper_bin[idx][0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.linalg.norm(p1-p2) < 1e-6:
					continue

				elif np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up, p1))
						sol_bin_down = np.vstack((sol_bin_down,p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up, p2 ))
						sol_bin_down = np.vstack((sol_bin_down,p1))
					break

	sol_upper              = np.vstack((sol_upper, sol_bin_up  ))
	sol_lower              = np.vstack((sol_lower, sol_bin_down))   

# sort the solutions
#########################
	direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_upper            = np.arccos (np.dot (direction, central_axis))
	sorted_theta_upper_idx = np.argsort (theta_upper)

	sol_upper = sol_upper [sorted_theta_upper_idx]
	sol_lower = sol_lower [sorted_theta_upper_idx]

	return sol_upper, sol_lower

######################################

def root_finder_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string, iterr, chi_ps, chi_pc, chi_sc):

	if scale_string == "phi_a":
		scale_a = 1e+6 * (10**(iterr))
		scale_b = 1
	elif scale_string == "phi_b":
		scale_a = 1
		scale_b = 1e+6 * (10**(iterr))
	else:
		print (f"Bad string provided: {scale_string}.")
		exit()

	direction     = (sol_lower[:,0:2] - center)/np.linalg.norm(sol_lower[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_lower   = np.arccos(np.dot(direction, central_axis))
	theta1        = theta_lower [max_ind]
	theta2        = theta_lower [max_ind+1]

	unsolved_lower_bin = binodal_lower[bad_idx][:, 0:2]
	direction_binodal  = (unsolved_lower_bin[:,0:2]-center)/np.linalg.norm(unsolved_lower_bin[:,0:2]-center, axis=1)[:, np.newaxis]
	theta_binodal      = np.arccos(np.dot(direction_binodal, central_axis))
	to_keep            = (theta_binodal > theta1) & (theta_binodal < theta2)

	unsolved_lower_bin = unsolved_lower_bin[to_keep]

	sol_bin_up    = np.empty((0,3))
	sol_bin_down  = np.empty((0,3))

	upper_guesses = np.linspace (sol_upper[max_ind], sol_upper[max_ind+1], 100)

	for idx in range (0, len(unsolved_lower_bin), 100):
		print (f"@ idx = {idx}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(phi[0], unsolved_lower_bin[idx][1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], unsolved_lower_bin[idx][1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], unsolved_lower_bin[idx][1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for iidx in range(len(upper_guesses)):
			root = fsolve (mu_equations, [unsolved_lower_bin[idx][0], upper_guesses[iidx][0]/scale_a, upper_guesses[iidx][1]/scale_b])

			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue

			else:
				p1 = np.array([root[0], unsolved_lower_bin[idx][1], 1-root[0]-unsolved_lower_bin[idx][1]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.linalg.norm(p1-p2) < 1e-6:
					continue

				elif np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up,  p1))
						sol_bin_down = np.vstack((sol_bin_down,p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up,  p2))
						sol_bin_down = np.vstack((sol_bin_down,p1))
					break

	for idx in range (0, len(unsolved_lower_bin), 100):
		print (f"@ idx = {idx}...", flush=True)
		def mu_equations (phi):
			eq1 = mu_a(unsolved_lower_bin[idx][0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(unsolved_lower_bin[idx][0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(unsolved_lower_bin[idx][0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for iidx in range(len(upper_guesses)):
			root = fsolve (mu_equations, [unsolved_lower_bin[idx][1], upper_guesses[iidx][0]/scale_a, upper_guesses[iidx][1]/scale_b])

			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue

			else:
				p1 = np.array([unsolved_lower_bin[idx][0], root[0], 1-root[0]-unsolved_lower_bin[idx][0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.linalg.norm(p1-p2) < 1e-6:
					continue

				elif np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)

					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up,  p1))
						sol_bin_down = np.vstack((sol_bin_down,p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up,  p2))
						sol_bin_down = np.vstack((sol_bin_down,p1))
					break

	sol_upper              = np.vstack((sol_upper, sol_bin_up  ))
	sol_lower              = np.vstack((sol_lower, sol_bin_down))

# sort the solutions
#########################
	direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_upper            = np.arccos  (np.dot (direction, central_axis))
	sorted_theta_upper_idx = np.argsort (theta_upper)

	sol_upper = sol_upper [sorted_theta_upper_idx]
	sol_lower = sol_lower [sorted_theta_upper_idx]

	return sol_upper, sol_lower

######################################

def binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc, vs, vp, vc, crit, center):

	try:
		df = pd.read_csv (dumpfile, sep='\s+', engine="python", skiprows=1, names=["dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2"])
		df = df.loc[df["dmu"]<5]
	except FileNotFoundError:
		print (f"File called {dumpfile} was not found. This was likely because pskelbin.py could not find reasonable guesses. Please check your parameters and inputs and try again.", flush=True)
		exit ()

	# def stab_crit (p_s, p_p, c_ps, c_pc, c_sc):
	# 	return (1/(vp*p_p) + 1/(vc*(1-p_s - p_p)) - 2 * c_pc) * (1/(vs*p_s) + 1/(vc*(1-p_s - p_p)) - 2 * c_sc) - (1/(vc*(1-p_s-p_p)) + c_ps - c_pc - c_sc) ** 2

# original point and the guess for the root... 
	phi_a_upper = df["phi_a1"].values; phi_a_lower = df["phi_a2"].values
	phi_b_upper = df["phi_b1"].values; phi_b_lower = df["phi_b2"].values
	phi_c_upper = df["phi_c1"].values; phi_c_lower = df["phi_c2"].values

	binodal_upper  = np.zeros ((phi_a_upper.shape[0],3))
	binodal_lower  = np.zeros ((phi_a_upper.shape[0],3))

	sol_upper = np.empty((0,3))
	sol_lower = np.empty((0,3))

	bad_idx  = []
	good_idx = []

# start partitioning along a certain axis
	central_axis   = (crit-center)/np.linalg.norm (crit-center)

	print (f"center = {center}", flush=True)
	print (f"central axis = {central_axis}", flush=True)

# f = open (args.boundary, 'w')
	print ("Start processing the dumpfile and find roots.", flush=True)
	print (f"I will be processing an array of size = {phi_a_upper.size}.", flush=True)

	for idx in range (len(phi_a_upper)):

		def mu_equations (phi):
			eq1 = mu_a(phi[0], phi_b_upper[idx]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], phi_b_upper[idx]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], phi_b_upper[idx]) - mu_c(phi[1], phi[2])

			return [eq1, eq2, eq3]

		root = fsolve (mu_equations, [phi_a_upper[idx], phi_a_lower[idx], phi_b_lower[idx]])
		binodal_upper [idx, :] = np.array([phi_a_upper[idx], phi_b_upper[idx], 1-phi_a_upper[idx]-phi_b_upper[idx]])
		binodal_lower [idx, :] = np.array([phi_a_lower[idx], phi_b_lower[idx], 1-phi_a_lower[idx]-phi_b_lower[idx]])

		# if the roots are "bad" roots, just write them out as bad
		if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
			bad_idx.append (idx)
			continue

		else:
			fa = [root[0], root[1]]
			fb = [phi_b_upper[idx], root[2]]
			fc = [1-root[0]-phi_b_upper[idx], 1-root[1]-root[2]]
			p1 = np.array([root[0], phi_b_upper[idx], 1-root[0]-phi_b_upper[idx]])
			p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

			if np.linalg.norm(p1-p2) < 1e-6:
				bad_idx.append (idx)

			elif args.ous and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
				bad_idx.append(idx)

			elif np.isnan(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				bad_idx.append (idx)

			elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				bad_idx.append (idx)

			elif np.sign(np.cross(central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
				bad_idx.append (idx)

			else:
				good_idx.append (idx)
				if np.cross (central_axis, p1[0:2]-center) >= 0:
					sol_upper = np.vstack ((sol_upper,p1))
					sol_lower = np.vstack ((sol_lower,p2))
				else:
					sol_upper = np.vstack ((sol_upper,p2))
					sol_lower = np.vstack ((sol_lower,p1))

	print (f"Length of bad_idx = {len(bad_idx)}")
	print (f"Length of good_idx = {len(good_idx)}")

	##################
	##################

	direction              = (sol_upper[:,0:2] - center) / np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
	theta_upper            = np.arccos (np.dot (direction, central_axis))
	sorted_theta_upper_idx = np.argsort (theta_upper)
	theta_upper            = theta_upper[sorted_theta_upper_idx]
	sol_upper              = sol_upper[sorted_theta_upper_idx]
	sol_lower              = sol_lower[sorted_theta_upper_idx]

# I have sorted the solutions 
# now, if the solution curve is sufficiently close, no need to perform more detailed searches -- so find maximum distances between points on the solution curves
# find differences between lower and upper curves
	diff_up       = np.linalg.norm(sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
	max_diff_up   = np.max(diff_up)

	diff_down     = np.linalg.norm(sol_lower[1:][:,0:2] - sol_lower[:-1][:,0:2], axis=1)
	max_diff_down = np.max(diff_down)

	max_diff_count = 0
	while (max_diff_up > 0.1 or max_diff_down > 0.1) and max_diff_count != args.dc:
		
		# if max_diff_up > max_diff_down:
		print (f"Running a search on the top half...", flush=True)
		# start running a finer search 
		difference    = np.linalg.norm (sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
		max_ind       = np.argmax (difference)

		#upper edges 
		u1 = sol_upper[max_ind][0:2]
		u2 = sol_upper[max_ind+1][0:2]
		umin = np.min([u1,u2])
		l1 = sol_lower[max_ind][0:2]
		l2 = sol_lower[max_ind+1][0:2]
		udist = abs(np.min(u1)-np.min(u2))
		ldist = abs(np.min(l1)-np.min(l2))

		print (f"u1 = {u1}, u2 = {u2}, l1 = {l1}, l2 = {l2}")
		print (f"udist = {udist}, ldist = {ldist}")
		if udist > ldist:
			print (f"In the scaling lower guesses zone...", flush=True)
			# scale lower guesses
			if   ( abs(l1[0] - l2[0]) < abs(l1[1]-l2[1]) ):
				scale_string = "phi_a"
				print (f"scale_string = {scale_string}")
				sol_upper, sol_lower = root_finder_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string, max_diff_count, chi_ps, chi_pc, chi_sc)

			elif ( abs(l1[0] - l2[0]) > abs(l1[1]-l2[1]) ):
				scale_string = "phi_b"
				print (f"scale_string = {scale_string}")
				sol_upper, sol_lower = root_finder_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string,max_diff_count, chi_ps, chi_pc, chi_sc)

		elif udist < ldist:
			print (f"In the scaling upper guesses zone...", flush=True)
			# scale upper guesses
			if   ( abs(u1[0] - u2[0]) < abs(u1[1]-u2[1]) ):
				scale_string = "phi_a"
				print (f"scale_string = {scale_string}")
				sol_upper, sol_lower = root_finder_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string, max_diff_count, chi_ps, chi_pc, chi_sc)

			elif ( abs(u1[0] - u2[0]) > abs(u1[1]-u2[1]) ):
				scale_string = "phi_b"
				print (f"scale_string = {scale_string}")
				sol_upper, sol_lower = root_finder_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string, max_diff_count, chi_ps, chi_pc, chi_sc)


		diff_up       = np.linalg.norm(sol_upper[1:] - sol_upper[:-1], axis=1)
		max_diff_up   = np.max(diff_up)

		diff_down     = np.linalg.norm(sol_lower[1:] - sol_lower[:-1], axis=1)
		max_diff_down = np.max(diff_down)
		max_diff_count += 1

		direction              = (sol_upper[:,0:2] - center) / np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
		theta_upper            = np.arccos (np.dot (direction, central_axis))
		sorted_theta_upper_idx = np.argsort (theta_upper)
		theta_upper            = theta_upper[sorted_theta_upper_idx]
		sol_upper              = sol_upper[sorted_theta_upper_idx]
		sol_lower              = sol_lower[sorted_theta_upper_idx]
		print (f"max_diff_count = {max_diff_count}.")
		# if max_diff_count == args.dc:
		#    break

##################################################

	print ("Broke out of initial solver. Time to refine...", flush=True)
# WE NOW HAVE A PRETTY CLEANED OUT BINODAL. 
# ALL THAT REMAINS IS TO REFINE IT TO MAKE IT LOOK NICE.

	diff_up       = np.linalg.norm(sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
	max_diff_up   = np.max(diff_up)

	diff_down     = np.linalg.norm(sol_lower[1:][:,0:2] - sol_lower[:-1][:,0:2], axis=1)
	max_diff_down = np.max(diff_down)

	print ("Being refining!", flush=True)
	max_refine_count = 0

	if args.rbs == -1:
		while_condition = (max_diff_up > 0.01 or max_diff_down > 0.01) and max_refine_count != args.rc
	else:
		while_condition = (len(sol_upper) < args.rbs) or ((max_refine_count != args.rc) and (max_diff_up > 0.01 or max_diff_down > 0.01))

	while while_condition: # (max_diff_up > 0.01 or max_diff_down > 0.01) and max_refine_count != args.rc:

		print (f"@ max_diff_up = {max_diff_up}, max_diff_down = {max_diff_down}...")
		nadded_rows = args.nar

		sol_upper, sol_lower   = refined_binodal_v4 (sol_upper, sol_lower, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc)

		direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
		theta_upper            = np.arccos (np.dot (direction, central_axis))
		sorted_theta_upper_idx = np.argsort (theta_upper)
		theta_upper            = theta_upper[sorted_theta_upper_idx]
		sol_upper              = sol_upper[sorted_theta_upper_idx]
		sol_lower              = sol_lower[sorted_theta_upper_idx]

		sol_upper, sol_lower   = refined_binodal_v5 (sol_upper, sol_lower, central_axis, nadded_rows, chi_ps, chi_pc, chi_sc)

		direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
		theta_upper            = np.arccos (np.dot (direction, central_axis))
		sorted_theta_upper_idx = np.argsort (theta_upper)
		theta_upper            = theta_upper[sorted_theta_upper_idx]
		sol_upper              = sol_upper[sorted_theta_upper_idx]
		sol_lower              = sol_lower[sorted_theta_upper_idx]

		diff_up       = np.linalg.norm(sol_upper[1:] - sol_upper[:-1], axis=1)
		max_diff_up   = np.max(diff_up)

		diff_down     = np.linalg.norm(sol_lower[1:] - sol_lower[:-1], axis=1)
		max_diff_down = np.max(diff_down)
		max_refine_count += 1
		print (f"max_refine_count = {max_refine_count}.")

		if args.rbs == -1:
			while_condition = (max_diff_up > 0.01 or max_diff_down > 0.01) and max_refine_count != args.rc
		else:
			while_condition = (len(sol_upper) < args.rbs) # and ((max_refine_count != args.rc) or (max_diff_up > 0.01 or max_diff_down > 0.01))
			print (f"len(sol_upper) = {len(sol_upper)}, args.rbs = {args.rbs}, (len(sol_upper) < args.rbs) = {len(sol_upper) < args.rbs}.", flush=True)
			print (f"while_condition = {while_condition}", flush=True)
#####

	ref_bin = [sol_upper, sol_lower]
	print ("This particular crit points should be well-populated.", flush=True)

# this is the binodal
	if args.ternary:
		ax.scatter (ref_bin[0][:,0], 1-ref_bin[0][:,0]-ref_bin[0][:,1], ref_bin[0][:,1], c='silver',     s=0.125, zorder=11)
		ax.scatter (ref_bin[1][:,0], 1-ref_bin[1][:,0]-ref_bin[1][:,1], ref_bin[1][:,1], c='darkred',    s=0.125, zorder=11)

	else:
		ax.scatter (ref_bin[0][:,0], ref_bin[0][:,1], c='silver',     s=0.125, zorder=11)
		ax.scatter (ref_bin[1][:,0], ref_bin[1][:,1], c='darkred',    s=0.125, zorder=11)
		ax.scatter (center[0], center[1], c='darkgreen', s=1, zorder=12)

	if args.tl:
		for i in range (0, len(ref_bin[0]), len(ref_bin[0])//args.td):
			if args.ternary:
				ax.plot    ([ref_bin[0][i,0],ref_bin[1][i,0]], \
							[1-ref_bin[0][i,0]-ref_bin[0][i,1], 1-ref_bin[1][i,0]-ref_bin[1][i,1]], \
							[ref_bin[0][i,1],ref_bin[1][i,1]], \
							lw=0.5, ls='--', markersize=0, zorder=10, c='skyblue')
			else:
				ax.plot    ([ref_bin[0][i,0],ref_bin[1][i,0]], \
						   [ref_bin[0][i,1],ref_bin[1][i,1]], \
						   lw=0.5, ls='--', markersize=0, zorder=10, c='skyblue')

	ff = open (boundaryfile, 'a')
	for i in range (len(ref_bin[0])):
		ff.write (f"{ref_bin[0][i][0]}|{ref_bin[0][i][1]}|{ref_bin[0][i][2]}|{ref_bin[1][i][0]}|{ref_bin[1][i][1]}|{ref_bin[1][i][2]}\n")

	return

############################
############################


if __name__=="__main__":

	start = time.time()

	###########################
	vs       = args.vs; va = args.vs
	vp       = args.vp; vb = args.vp
	vc       = args.vc; vc = args.vc
	chi_ps   = args.chi_ps
	chi_pc   = args.chi_pc
	chi_sc   = args.chi_sc
	dumpfile = args.dumpfile
	############################
	lsize = 3
	font = {'color':  'black',
		'weight': 'normal',
		'size': lsize}

	fig = plt.figure(num=1, figsize=(8,8))
	if args.ternary:
		ax = fig.add_subplot (projection="ternary")
	else:
		ax  = plt.axes ()
	###################

	# solution in terms of phi_s
	discriminant_s = lambda phi_s: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc)) + (vp - 2*phi_s*vp *vs *chi_ps + vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc) ) ) )**2 

	denom_s    = lambda phi_s: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))
	prefac_s   = lambda phi_s: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp * (2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc) ) )
	root_up_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) + np.sqrt(discriminant_s(phi_s)))
	root_lo_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) - np.sqrt(discriminant_s(phi_s)))

    # solution in terms of phi_p
	discriminant_p = lambda phi_p: -4*vc*vs*(phi_p*vp+(-1+phi_p)*vc*(-1+2*phi_p*vp*chi_pc))*(2*chi_sc+phi_p*vp*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc)))+(vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2-2*chi_pc*(1+vs*(chi_ps+chi_sc))))))**2
	denom_p        = lambda phi_p: 1/(-2*vc*vs*(2*chi_sc+phi_p*vp*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))))
	prefac_p       = lambda phi_p: vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2-2*chi_pc*(1+vs*(chi_ps+chi_sc)))))
	root_up_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)+np.sqrt(discriminant_p(phi_p)))
	root_lo_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)-np.sqrt(discriminant_p(phi_p)))
    ###################

	roots_up, roots_down = ternary.find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s)
	crits = np.vstack ((roots_up, roots_down))
	crits = ternary.remove_close_rows (crits, 1e-3)

	# plotting out the critical points
	if args.ternary:
		ax.scatter (crits[:,0], 1-crits[:,0]-crits[:,1], crits[:,1], color='k', edgecolors='greenyellow', s=1, zorder=15)
	else:
		ax.scatter (crits[:,0], crits[:,1], color='k', edgecolors='greenyellow', s=1)

	if len(crits) == 2:
		pass
	elif len(crits) == 0:
		print(f"There are no critical points. Exiting...")
		exit ()
	elif len(crits) == 3:
		print(f"Number of critical points is {len(crits)} == 3. Currently, we do not have the ability to solve for such diagrams. Here we go...")
		pass
	elif len(crits) > 4:
		print(f"Number of critical points is {len(crits)} > 3. Currently, we do not have the ability to solve for such diagrams. Here we go...")
		pass

	print ("Number of critical points = ",len(crits), flush=True)
	print ("critical points = \n",crits, flush=True)
	print ("Begin creating the meshes and painting the ternary diagram...", flush=True)
	p_s_space = np.arange (0.001, 1-0.001, 0.001)
	p_s = np.repeat (p_s_space, len(p_s_space))

	p_p = np.zeros (p_s.shape)
	for i in range (len(p_s_space)):
		p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

	vals = ternary.stab_crit (p_s, p_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)

	to_keep = ~np.isnan(vals)

	vals = vals[to_keep]
	p_s  = p_s [to_keep]
	p_p  = p_p [to_keep]

	vmax = np.max (vals)
	vmin = np.min (vals)

	if np.sign (vmax) == np.sign (vmin):
		if np.sign (vmax) >=0:
			vmin = -vmax
		else:
			vmax = -vmin

	print (f"vmin = {vmin}, vmax = {vmax}", flush=True)

	if len(vals) == 0:
		print(f"There is no spinodal region.")
		exit()

	norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
	cols = cm.bwr (norm (vals))

	# Plot the points
	p_c = 1 - p_s - p_p
	if args.ternary:
		ax.scatter (p_s, p_c, p_p, s=1, color=cols, zorder=0, clip_on=True)
	else:
		ax.scatter  (p_s, p_p, s=0.01, color=cols, zorder=0, clip_on=True)

	print ("Painted the ternary diagram!", flush=True)

	print ("We have plotted the spinodal region!\n", flush=True)
	print ("###########################################################\n", flush=True)
	print ("Start binodal plotting...\n", flush=True)

	mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_sc + phi_b * (1-phi_a-phi_b) * (chi_ps + chi_sc - chi_pc) ) 
	mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_pc + phi_a * (1-phi_a-phi_b) * (chi_ps + chi_pc - chi_sc) )
	mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )

	ff = open (boundaryfile, 'w')
	ff.write  ("phi_s_top|phi_p_top|phi_c_top|phi_s_bot|phi_p_bot|phi_c_bot\n")
	ff.close  ()
	normals  = []
	tangents = []
	for idx,crit in enumerate(crits):
		print (f"@ crit point = {crit}...")
		tangent_to_crit = tangent.tangent2  (vs, vc, vp, crit[0], crit[1], chi_pc, chi_ps, chi_sc, root_up_s, root_lo_s)
		normal_slope    = -1/tangent_to_crit
		tangents.append (tangent_to_crit)
		if len(crits) == 1:
			center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit
			binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc, va, vb, vc, crit, center)
		elif len(crits) == 2:
			if args.normal:
				center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit
				binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc, va, vb, vc, crit, center)
				if idx == 1:
					pseudo_center = np.mean (crits, axis=0)
					direction = np.mean(tangents, axis=0)/np.linalg.norm(np.mean(tangents, axis=0))
					pseudo_crit = center + direction
					print ("The pseudo_center is {pseudo_center}, the pseudocrit point is {pseudo_crit}.", flush=True)
					binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc,va, vb, vc, pseudo_crit, pseudo_center)
			else:
				center  = np.mean(crits, axis=0)
				binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc, va, vb, vc, crit, center)
				break
		else:
			center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit
			binodal_plotter (fig, ax, dumpfile, chi_ps, chi_pc, chi_sc, va, vb, vc, crit, center)

	# do a final sort
	print ("All crit points have been addressed.")
	print ("Done with binodal plotting!", flush=True)

	ternary.embelish(ax, args.ternary)
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
		plt.savefig (f"binodal_tern-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)

	else:
		plt.savefig (f"binodal_reg-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)

	stop = time.time()
	print (f"Time for computation is {stop-start} seconds.", flush=True)

