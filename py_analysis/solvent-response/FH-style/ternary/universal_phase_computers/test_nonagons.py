import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import tangent
import mpltern
import pickle
import phase
from shapely.geometry import LineString, MultiLineString, Point

import argparse
parser = argparse.ArgumentParser(description='Create a skeleton solution for the binodal. This is a memory-intensive computation.'         )
parser.add_argument('--chisc',               metavar='chi_sc',  dest='chi_sc',        type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',               metavar='chi_ps',  dest='chi_ps',        type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',               metavar='chi_pc',  dest='chi_pc',        type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',                   metavar='vs',      dest='vs',            type=float, action='store', help='specific volume of solvent.'  )
parser.add_argument('-vc',                   metavar='vc',      dest='vc',            type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',                   metavar='vp',      dest='vp',            type=float, action='store', help='specific volume of polymer.'  )
parser.add_argument('--final-binodal',       metavar='FB',      dest='fb',            type=str,   action='store', help='name of pickle file to dump the BINODAL object in.')
parser.add_argument('--search-density',      metavar='SD',      dest='sd',            type=int,   action='store', help='density of points sampled for stability plot (default: 500).',                          default=500 )
parser.add_argument('--island-stable-pkl',   metavar='SPKL',    dest='spkl',          type=str,   action='store', help='extract information about the stable islands from the pickle file (default: None).',    default=None)
parser.add_argument('--island-unstable-pkl', metavar='UPKL',    dest='upkl',          type=str,   action='store', help='extract information about the unstable islands from the pickle file (default: None).',  default=None)
parser.add_argument('--crit-pkl',            metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                                default=None)
parser.add_argument('--binodal-pkl',         metavar='bpkl',    dest='bpkl',          type=str,   action='store', help='enter name of file with all the information about the binodals (default: None).',       default=None)
parser.add_argument('--plot-edges',     dest='pe',         action='store_true',  help='plot the edges of the spinodal.',     default=False)
parser.add_argument('--plot-crits',     dest='pc',         action='store_true',  help='plot the critical points.'      ,     default=False)
parser.add_argument('--no-rtw',         dest='nrtw',       action='store_true',  help="Dont print out the runtime warning.", default=False)
parser.add_argument('--img',            dest='img',        action='store',       type=str, default="None", help="name of image to be created.")
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

def transform_islands(islands):
	hull_paths = []
	for idx, island in enumerate(islands):
		islands[idx] = np.array(np.vstack([0.001+(0.999-0.001)/args.sd*island[:,1], 0.001+(0.999-0.001)/args.sd*island[:,0]])).T 
		hull = ConvexHull(islands[idx])
		hull_paths.append(Path(islands[idx][hull.vertices]))

	return hull_paths

#########################################

def find_central_island(islands):
	c_dists = []
	center  = np.array([1/3,1/3])
	for i in range(len(islands)):
		island_c = np.mean(islands[i], axis=0)
		c_dists.append(np.linalg.norm(center-island_c))
	
	c_idx = np.argmin(c_dists)

	return c_idx 

#########################################

def truncate_curve_at_point(curve, point):
	"""
	Truncate a curve represented by a NumPy array at a specified point.

	Parameters:
	- curve (numpy.ndarray): Array representing the curve.
	- point (numpy.ndarray): Coordinates of the point to truncate the curve at.

	Returns:
	- truncated_curve (numpy.ndarray): Truncated curve up to and including the specified point.
	"""

	# Calculate distances from each point in curve to point
	distances_P = np.linalg.norm(curve - point, axis=1)

	# Find the index where the distances change signs
	min_index = np.argmin(distances_P)

	# check which side of min_index this point is
	if np.linalg.norm(curve[min_index+1] - point) > np.linalg.norm(curve[min_index-1] - point):
		insertion_index = min_index
		# print(f"curve inserts = {curve[min_index-1], point, curve[min_index]}")
	else:
		insertion_index = min_index + 1
		# print(f"curve inserts = {curve[min_index], point, curve[min_index+1]}")

	return insertion_index

#########################################

def clean_and_sort(sol1, sol2, crit, norm_vec):

	adj_sol1 = (sol1[:,0:2]-crit[0:2])/np.linalg.norm(sol1[:, 0:2]-crit[0:2], axis=1).reshape(-1,1)
	adj_sol2 = (sol2[:,0:2]-crit[0:2])/np.linalg.norm(sol2[:, 0:2]-crit[0:2], axis=1).reshape(-1,1)
	signs1   = np.sign(np.cross(norm_vec, adj_sol1))
	signs2   = np.sign(np.cross(norm_vec, adj_sol2))

	# find indices where signs are +1 and -1
	good_indices_1 = np.where((signs1== 1) & (signs2==-1))[0]
	good_indices_2 = np.where((signs1==-1) & (signs2== 1))[0]

	pos_sol = np.vstack((sol1[good_indices_1], sol2[good_indices_2]))
	neg_sol = np.vstack((sol2[good_indices_1], sol1[good_indices_2]))
	
	dists   = np.linalg.norm(pos_sol[:,0:2]-crit[0:2], axis=1)
	pos_sol = pos_sol[np.argsort(dists)]
	neg_sol = neg_sol[np.argsort(dists)]

	return neg_sol, pos_sol

#########################################

def check_sidedness(P, arm, crit):
	adj_arm    = (arm[:, 0:2] - crit[0:2])/np.linalg.norm(arm[:, 0:2] - crit[0:2], axis=1).reshape(-1,1)
	tang_slope = tangent.tangent2(P.vs, P.vc, P.vp, crit[0], crit[1], P.chi_pc, P.chi_ps, P.chi_sc, P.spinodal.root_up_s, P.spinodal.root_lo_s)
	norm_slope = -1/tang_slope 
	norm_vec   = np.array([1, norm_slope])/np.sqrt(1+norm_slope**2)
	signs      = np.sign(np.cross(norm_vec, adj_arm))
	if len(np.unique(signs)) == 1:
		return 0
	else:
		return 1

#########################################

def check_sidedness_list(P, arm, critlist):
	bools = []
	for crit in critlist:
		bools.append(check_sidedness(P, arm, crit))
	
	if len(np.unique(bools)) == 1:
		return True
	else:
		return False

#########################################
	
def find_overlapping_sides(P, arm, critlist):
	signs_vecs = []
	for crit in critlist:
		adj_arm    = (arm[:, 0:2] - crit[0:2])/np.linalg.norm(arm[:, 0:2] - crit[0:2], axis=1).reshape(-1,1)
		tang_slope = tangent.tangent2(P.vs, P.vc, P.vp, crit[0], crit[1], P.chi_pc, P.chi_ps, P.chi_sc, P.spinodal.root_up_s, P.spinodal.root_lo_s)
		norm_slope = -1/tang_slope 
		norm_vec   = np.array([1, norm_slope])/np.sqrt(1+norm_slope**2)
		signs      = np.sign(np.cross(norm_vec, adj_arm))
		signs_vecs.append(signs) 

	# now, look where signs are not the same for two arrays 
	to_extract = []
	for i in range(len(signs_vecs)):
		for j in range(i+1, len(signs_vecs)):
			unequal = (signs_vecs[i] != signs_vecs[j]) 
			to_extract.append(unequal)

	return

#########################################

def remove_intersecting_segments(sol1, sol2):

	# Convert the line segments to Shapely LineString objects
	lines = np.array([LineString([point1, point2]) for point1, point2 in zip(sol1[:,0:2], sol2[:,0:2])])

	# Create a matrix of intersections using broadcasting
	intersects_matrix = np.triu(np.array([[lines[i].intersects(lines[j]) for j in range(len(lines))] for i in range(len(lines))]), k=1)

	# Find indices of non-intersecting line segments
	non_intersecting_indices = ~np.any(intersects_matrix, axis=0)

	# Extract non-intersecting line segments
	non_intersecting_sol1 = sol1[non_intersecting_indices]
	non_intersecting_sol2 = sol2[non_intersecting_indices]

	return non_intersecting_sol1, non_intersecting_sol2

#########################################

def add_rows(array, M, idx):
	insert_rows = np.linspace(array[idx], array[idx+1], M)
	# pop_array   = np.insert(array, idx+1, insert_rows[1:-1], axis=0)
	# return pop_array
	return insert_rows

#########################################

def solve_within(arm_neg, arm_pos, P, center, central_axis):

	b_neg = np.empty((0,3))
	b_pos = np.empty((0,3))

	for idx, pt in enumerate(arm_neg):
		def mu_equations(phi):
			eq1 = P.sym_mu_ps.delta_mu_s(pt[0], phi[0], phi[1], phi[2]) 
			eq2 = P.sym_mu_ps.delta_mu_p(pt[0], phi[0], phi[1], phi[2]) 
			eq3 = P.sym_mu_ps.delta_mu_c(pt[0], phi[0], phi[1], phi[2]) 
			return [eq1, eq2, eq3]

		for tidx, tpt in enumerate(arm_pos):
			root = fsolve(mu_equations, [pt[1], tpt[0], tpt[1]])
			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue
			else:
				p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])
				# print(f"p1 = {p1}, p2 = {p2}")
				if np.isnan(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					# print(f"Nan points.")
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					# print(f"Inf points.")
					continue

				elif ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) < 0 or ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) < 0:
					# print(f"Unstable points.")
					continue 

				elif np.linalg.norm(p1[0:2]-p2[0:2]) < 1e-6:
					# print(f"Too close.")
					continue

				else:
					if np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) == np.sign(np.cross(central_axis, p2[0:2]-center[0:2])):
						# print(f"Same side.")
						continue
					elif np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) >=0 and np.sign(np.cross(central_axis, p2[0:2] - center[0:2])) <= 0:
						# print(f"Good!.")
						b_neg = np.vstack((b_neg, p2))
						b_pos = np.vstack((b_pos, p1))
					elif np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) <=0 and np.sign(np.cross(central_axis, p2[0:2] - center[0:2])) >= 0:
						# print(f"Good!")
						b_neg = np.vstack((b_neg, p1))
						b_pos = np.vstack((b_pos, p2))
					else:
						# print(f"Something else...")
						continue
					break

	return b_neg, b_pos
				
#########################################

def add_and_solve(arm_neg, arm_pos, P, center, central_axis, M):

	# find point with greatest gap
	print(f"arm_neg.shape = {arm_neg.shape}", flush=True)
	print(f"arm_pos.shape = {arm_pos.shape}", flush=True)
	dist_neg     = np.linalg.norm(np.diff(arm_neg, axis=0), axis=1)
	max_dist_neg = np.argmax(dist_neg)
	dist_pos     = np.linalg.norm(np.diff(arm_pos, axis=0), axis=1)
	max_dist_pos = np.argmax(dist_pos)

	# print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
	# print(f"Negative points are: {arm_neg[max_dist_neg], arm_neg[max_dist_neg+1]}")
	# print(f"Positive points are: {arm_pos[max_dist_pos], arm_pos[max_dist_pos+1]}")
	iterx = 0

	while dist_neg[max_dist_neg] > 0.0005 or dist_pos[max_dist_pos] > 0.0005:

		# print(f"Negative points are: {arm_neg[max_dist_neg], arm_neg[max_dist_neg+1]}")
		# print(f"Positive points are: {arm_pos[max_dist_pos], arm_pos[max_dist_pos+1]}")

		if dist_neg[max_dist_neg] > dist_pos[max_dist_pos] and dist_neg[max_dist_neg] > 0.0005:
			# print("Filling in for negative distances...", flush=True)
			new_arm_neg = add_rows(arm_neg, M, max_dist_neg)
			new_arm_pos = add_rows(arm_pos, M, max_dist_neg) 

			# print(f"the new arm neg is =\n{new_arm_neg[:5], new_arm_neg[-5:]} ")
			# print(f"the new arm pos is =\n{new_arm_pos[:5], new_arm_pos[-5:]} ")

			# b_neg, b_pos  = solve_within(new_arm_neg[max_dist_neg+1:max_dist_neg+1+M], new_arm_pos[max_dist_neg+1:max_dist_neg+1+M], P, center, central_axis)
			b_neg, b_pos  = solve_within(new_arm_neg, new_arm_pos, P, center, central_axis)

			# print(f"b_neg = \n{b_neg}, b_pos = \n{b_pos}")

			arm_neg = np.vstack((arm_neg, b_neg))
			arm_pos = np.vstack((arm_pos, b_pos))
			arm_neg, arm_pos = clean_and_sort(arm_neg, arm_pos, center, central_axis)


		elif dist_pos[max_dist_pos] > dist_neg[max_dist_neg] and dist_neg[max_dist_neg] > 0.0005:
			# print("Filling in for positive distances...", flush=True)
			new_arm_neg = add_rows(arm_neg, M, max_dist_pos)
			new_arm_pos = add_rows(arm_pos, M, max_dist_pos)

			# print(f"the new arm neg is =\n{new_arm_neg[:5], new_arm_neg[-5:]} ")
			# print(f"the new arm pos is =\n{new_arm_pos[:5], new_arm_pos[-5:]} ")

			# b1, b2  = solve_within(new_arm_neg[max_dist_pos+1:max_dist_pos+1+M], new_arm_pos[max_dist_pos+1:max_dist_pos+1+M], P, center, central_axis)
			b_neg, b_pos  = solve_within(new_arm_neg, new_arm_pos, P, center, central_axis)
			arm_neg = np.vstack((arm_neg, b_neg))
			arm_pos = np.vstack((arm_pos, b_pos))
			arm_neg, arm_pos = clean_and_sort(arm_neg, arm_pos, center, central_axis)

		# find point with greatest gap
		arm_neg, keep = ternary.remove_close_rows(arm_neg, 1e-12)
		arm_pos       = arm_pos[keep]

		dist_neg      = np.linalg.norm(np.diff(arm_neg, axis=0), axis=1)
		max_dist_neg  = np.argmax(dist_neg)
		dist_pos      = np.linalg.norm(np.diff(arm_pos, axis=0), axis=1)
		max_dist_pos  = np.argmax(dist_pos)


		print(f"iterx = {iterx}", flush=True)
		# print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
		iterx += 1
		if iterx > 5:
			break

	return arm_neg, arm_pos 

#########################################

def distance_from_axis(point, axial_point, axis):
	r     = np.linalg.norm(point[0:2]-axial_point[0:2])
	rhat  = (point[0:2] - axial_point[0:2])/np.linalg.norm(point[0:2] - axial_point[0:2])
	theta = np.arccos(np.sum(rhat[0:2]*axis[0:2]))
	distance = r*np.sin(theta)
	return distance 

##########################################

def solve_central(arm_neg, arm_pos, center, central_axis):

	b_neg = np.empty((0,3))
	b_pos = np.empty((0,3))

	for idx, pt in enumerate(arm_neg):
		def mu_equations(phi):
			eq1 = P.sym_mu_ps.delta_mu_s(pt[0], phi[0], phi[1], phi[2]) 
			eq2 = P.sym_mu_ps.delta_mu_p(pt[0], phi[0], phi[1], phi[2]) 
			eq3 = P.sym_mu_ps.delta_mu_c(pt[0], phi[0], phi[1], phi[2]) 
			return [eq1, eq2, eq3]

		for tidx, tpt in enumerate(arm_pos):
			root = fsolve(mu_equations, [pt[1], tpt[0], tpt[1]])
			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue
			else:
				p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])
				# print(f"p1 = {p1}, p2 = {p2}")
				if np.isnan(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					# print(f"Nan points.")
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					# print(f"Inf points.")
					continue

				elif ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) < 0 or ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) < 0:
					# print(f"Unstable points.")
					continue 

				elif np.linalg.norm(p1[0:2]-p2[0:2]) < 1e-6:
					# print(f"Too close.")
					continue

				else:
					if np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) == np.sign(np.cross(central_axis, p2[0:2]-center[0:2])):
						# print(f"Same side.")
						continue
					elif np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) >=0 and np.sign(np.cross(central_axis, p2[0:2] - center[0:2])) <= 0:
						# print(f"Good!")
						b_neg = np.vstack((b_neg, p2))
						b_pos = np.vstack((b_pos, p1))
					elif np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) <=0 and np.sign(np.cross(central_axis, p2[0:2] - center[0:2])) >= 0:
						# print(f"Good!")
						b_neg = np.vstack((b_neg, p1))
						b_pos = np.vstack((b_pos, p2))
					else:
						# print(f"Something else...")
						continue
					break

	return b_neg, b_pos

##########################################

def populate_center(BINODALS, uidx, M=50):

	nv_c = BINODALS["crit_info"][BINODALS["groupings"][uidx]["center"  ]["idx"]]["norm_vec"]
	nv_p = BINODALS["crit_info"][BINODALS["groupings"][uidx]["positive"]["idx"]]["norm_vec"]
	nv_n = BINODALS["crit_info"][BINODALS["groupings"][uidx]["negative"]["idx"]]["norm_vec"]

	cp_c = P.crits[BINODALS["groupings"][uidx]["center"  ]["idx"]]
	cp_p = P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]] 
	cp_n = P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]] 

	Bneg     = BINODALS["groupings"][uidx]["negative"]["binodals"]
	Bpos     = BINODALS["groupings"][uidx]["positive"]["binodals"]
	Bc       = BINODALS["groupings"][uidx]["center"  ]["binodals"]


	d_neg = distance_from_axis(Bneg[1][-1], cp_c, nv_c)
	d_pos = distance_from_axis(Bpos[0][-1], cp_c, nv_c)
	
	iterx = 0
	while d_neg > 0.0005 and d_pos > 0.0005:
		
		# treating the negative side 
		to_insert_pos = np.linspace(Bneg[1][-1], Bpos[0][-1], M)
		to_insert_neg = np.linspace(Bneg[0][-1], Bc[0][0],    M)

		arm_neg, arm_pos = solve_central(to_insert_pos, to_insert_neg, cp_n, nv_n) 
		# only keep points to the negative of the central axis 
		adj_arm_pos = (cp_c[0:2]-arm_pos[:,0:2])/np.linalg.norm(cp_c[0:2]-arm_pos[:,0:2], axis=1).reshape(-1,1)
		signs       = np.sign(np.cross(nv_c, adj_arm_pos))
		arm_neg     = arm_neg[signs<=0]
		arm_pos     = arm_pos[signs<=0]
		BINODALS["groupings"][uidx]["negative"]["binodals"][0] = \
			np.vstack((BINODALS["groupings"][uidx]["negative"]["binodals"][0], arm_neg))
		BINODALS["groupings"][uidx]["negative"]["binodals"][1] = \
			np.vstack((BINODALS["groupings"][uidx]["negative"]["binodals"][1], arm_pos))
		BINODALS["groupings"][uidx]["negative"]["binodals"][0], BINODALS["groupings"][uidx]["negative"]["binodals"][1] = \
			clean_and_sort(BINODALS["groupings"][uidx]["negative"]["binodals"][0], BINODALS["groupings"][uidx]["negative"]["binodals"][1], cp_n, nv_n) 

		d_neg = distance_from_axis(BINODALS["groupings"][uidx]["negative"]["binodals"][1][-1], cp_c, nv_c)
		
		# moving on to the positive side 
		to_insert_neg = np.linspace(Bpos[0][-1], Bneg[1][-1], M)
		to_insert_pos = np.linspace(Bpos[1][-1], Bc[1][0],    M)

		arm_neg, arm_pos = solve_central(to_insert_pos, to_insert_neg, cp_p, nv_p) 
		# only keep points to the negative of the central axis 
		adj_arm_neg = (cp_c[0:2]-arm_neg[:,0:2])/np.linalg.norm(cp_c[0:2]-arm_neg[:,0:2], axis=1).reshape(-1,1)
		signs       = np.sign(np.cross(nv_c, adj_arm_neg))
		arm_neg     = arm_neg[signs>=0]
		arm_pos     = arm_pos[signs>=0]
		BINODALS["groupings"][uidx]["positive"]["binodals"][0] = \
			np.vstack((BINODALS["groupings"][uidx]["positive"]["binodals"][0], arm_neg))
		BINODALS["groupings"][uidx]["positive"]["binodals"][1] = \
			np.vstack((BINODALS["groupings"][uidx]["positive"]["binodals"][1], arm_pos))
		BINODALS["groupings"][uidx]["positive"]["binodals"][0], BINODALS["groupings"][uidx]["positive"]["binodals"][1] = \
			clean_and_sort(BINODALS["groupings"][uidx]["positive"]["binodals"][0], BINODALS["groupings"][uidx]["positive"]["binodals"][1], cp_p, nv_p) 

		d_pos = distance_from_axis(BINODALS["groupings"][uidx]["positive"]["binodals"][1][-1], cp_c, nv_c)
		iterx += 1
		if iterx > 10:
			print("Breaking out of the center populater...", flush=True)
			break

				
		'''
		elif probe == "positive":
			
			iterx = 0
			while d_pos > 0.0005:
				to_insert_neg = np.linspace(Bpos[0][-1], Bneg[1][-1], M)
				to_insert_pos = np.linspace(Bpos[1][-1], Bc[1][0],    M)

				arm_neg, arm_pos = solve_central(to_insert_pos, to_insert_neg, cp_p, nv_p) 
				# only keep points to the negative of the central axis 
				adj_arm_neg = (cp_c[0:2]-arm_neg[:,0:2])/np.linalg.norm(cp_c[0:2]-arm_neg[:,0:2], axis=1).reshape(-1,1)
				signs       = np.sign(np.cross(nv_c, adj_arm_neg))
				arm_neg     = arm_neg[signs>=0]
				arm_pos     = arm_pos[signs>=0]
				BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], arm_neg))
				BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], arm_pos))
				BINODALS["groupings"][uidx][probe]["binodals"][0], BINODALS["groupings"][uidx][probe]["binodals"][1] = \
					clean_and_sort(BINODALS["groupings"][uidx][probe]["binodals"][0], BINODALS["groupings"][uidx][probe]["binodals"][1], cp_p, nv_p) 

				d_pos = distance_from_axis(BINODALS["groupings"][uidx][probe]["binodals"][1][-1], cp_c, nv_c)
				iterx += 1
				if iterx > 10:
					break
				# break
		'''
	return


##########################################

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

	if args.critpkl == None:	
		P.spinodal.obtain_crits()
	else:
		try:
			f = open(args.critpkl, 'rb')
			crits = pickle.load(f)
			P.spinodal.crits = crits
			f.close()
		except:
			P.spinodal.obtain_crits()
			f = open(args.critpkl, 'wb')
			pickle.dump(P.spinodal.crits, f)
			f.close()
	
	P.crits = P.spinodal.crits
	print(f"crits = \n{P.crits}", flush=True)
	print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	# P.tangent_tracing(ax)
	# print("Plotted out the tangent trace!", flush=True)
	
	# extract island information
	f = open(args.spkl, 'rb')
	stable_islands = pickle.load(f)
	f.close()

	f = open(args.upkl, 'rb')
	unstable_islands = pickle.load(f)
	f.close()

	#=================================
	# create all the volume fractions
	p_s_space = np.arange(0.001, 1-0.001, 0.001)
	p_s = np.repeat(p_s_space, len(p_s_space))

	p_p = np.zeros(p_s.shape)
	for i in range (len(p_s_space)):
		p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))
	
	# curate all the volume fractions
	vals = ternary.stab_crit(p_s, p_p, P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)
	to_keep = ~np.isnan(vals) & ~np.isinf(vals) & (vals > 0)

	p_s = p_s[to_keep]
	p_p = p_p[to_keep]
	phi = np.array([p_s, p_p]).T
	# print(f"phi = {phi}", flush=True)
	#=================================

	#=================================
	unstable_centers = []
	
	hull_paths_u = transform_islands(unstable_islands)
	# print(f"ui = {unstable_islands}", flush=True)

	for uidx, ui in enumerate(unstable_islands):
		unstable_centers.append(np.mean(ui, axis=0))
	unstable_centers = np.array(unstable_centers)
	# print(f"unstable = {unstable_centers}", flush=True)
	ax.scatter(unstable_centers[:,0], 1-unstable_centers[:,0]-unstable_centers[:,1], unstable_centers[:,1], c='orange', s=2, zorder=200)
	#=================================

	print(f"Number of stable islands is {len(stable_islands)} and unstable islands is {len(unstable_islands)}...", flush=True)
	
	
	f = open('binodals_test.pkl', 'rb')
	BINODALS = pickle.load(f)
	f.close() 
	uidx=0

	# run the command
	# populate_center(BINODALS, uidx, 50)
	# probe = "negative"
	# amped_neg_sol, amped_pos_sol = add_and_solve(BINODALS["groupings"][uidx][probe]["binodals"][0], BINODALS["groupings"][uidx][probe]["binodals"][1], P, P.crits[BINODALS["groupings"][uidx][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][uidx][probe]["idx"] ]["norm_vec"], 50)
	# BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][uidx][probe]["idx"]], amped_neg_sol))
	# BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][uidx][probe]["idx"]], amped_pos_sol))

	# unite the positive and negative binodals
	b_neg = BINODALS["groupings"][uidx]["negative"]["binodals"]
	b_pos = BINODALS["groupings"][uidx]["positive"]["binodals"]
	b_c   = BINODALS["groupings"][uidx]["center"  ]["binodals"]

	ax.scatter(b_neg[0][:,0], b_neg[0][:,2], b_neg[0][:,1], c='black', s=0.5)
	ax.scatter(b_neg[1][:,0], b_neg[1][:,2], b_neg[1][:,1], c='white', s=0.5)
	
	ax.scatter(b_pos[0][:,0], b_pos[0][:,2], b_pos[0][:,1], c='black', s=0.5)
	ax.scatter(b_pos[1][:,0], b_pos[1][:,2], b_pos[1][:,1], c='white', s=0.5)

	ax.scatter(b_c[0][:,0], b_c[0][:,2], b_c[0][:,1], c='black', s=0.5)
	ax.scatter(b_c[1][:,0], b_c[1][:,2], b_c[1][:,1], c='white', s=0.5)

	fig.savefig(args.img, dpi=1200, bbox_inches='tight')
	exit()
	
	# f = open("nonagons.pkl", 'wb')
	# pickle.dump(BINODALS, f)
	# f.close()
	
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

