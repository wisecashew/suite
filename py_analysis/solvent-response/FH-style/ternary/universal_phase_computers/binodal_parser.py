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
import ternary
import tangent
import mpltern
import pickle
import phase
import blibrary
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

import argparse
parser = argparse.ArgumentParser(description='Create a skeleton solution for the binodal. This is a memory-intensive computation.'         )
parser.add_argument('--chisc',               metavar='chi_sc',  dest='chi_sc',        type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',               metavar='chi_ps',  dest='chi_ps',        type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',               metavar='chi_pc',  dest='chi_pc',        type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',                   metavar='vs',      dest='vs',            type=float, action='store', help='specific volume of solvent.'  )
parser.add_argument('-vc',                   metavar='vc',      dest='vc',            type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',                   metavar='vp',      dest='vp',            type=float, action='store', help='specific volume of polymer.'  )
parser.add_argument('--final-binodal',       metavar='FB',      dest='fb',            type=str,   action='store', help='name of pickle file to dump the BINODAL object in.')
parser.add_argument('--database',            metavar='DB',      dest='db',            type=str,   action='store', help='name of database for this system.')
parser.add_argument('--island-stable-pkl',   metavar='SPKL',    dest='spkl',          type=str,   action='store', help='extract information about the stable islands from the pickle file (default: None).',   default=None)
parser.add_argument('--island-unstable-pkl', metavar='UPKL',    dest='upkl',          type=str,   action='store', help='extract information about the unstable islands from the pickle file (default: None).', default=None)
parser.add_argument('--crit-pkl',            metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                                default=None)
parser.add_argument('--binodal-pkl',         metavar='bpkl',    dest='bpkl',          type=str,   action='store', help='enter name of file with all the information about the binodals (default: None).',       default=None)
parser.add_argument('--plot-edges',     dest='pe',         action='store_true',  help='plot the edges of the spinodal.',     default=False)
parser.add_argument('--plot-crits',     dest='pc',         action='store_true',  help='plot the critical points.'      ,     default=False)
parser.add_argument('--plot-binodals',  dest='pb',         action='store_true',  help='plot the binodal points.'       ,     default=False)
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

# def find_truncation_point(curve, point):

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
	pos_indices = np.where((signs1== 1) & (signs2==-1))[0]
	neg_indices = np.where((signs1==-1) & (signs2== 1))[0]

	# use indices to separate rows into new arrays
	pos_sol = np.vstack((sol1[pos_indices], sol2[neg_indices]))
	neg_sol = np.vstack((sol2[pos_indices], sol1[neg_indices]))

	dists = np.linalg.norm(pos_sol[:,0:2]-crit[0:2], axis=1)
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
	# pop_array       = np.insert(array, idx+1, insert_rows[1:-1], axis=0)
	# return pop_array
	return insert_rows

#########################################

def find_closest_line_segment(A, B, P):

	# Calculate the vector from A to P
	AP = P - A
	dist_AP = np.linalg.norm(AP, axis=1)

	# calculate the vector from A to B
	AB = B - A
	dist_AB = np.linalg.norm(AB, axis=1)

	# calculate the angle between AP and AB
	theta = np.arccos( (np.sum(AP*AB, axis=1))/(dist_AP * dist_AB) )

	# get the perpendicular distance
	PP_dash = dist_AP * np.sin(theta)

	# get the smallest distance
	closest = np.argmin(PP_dash)

	return closest

#########################################

def identify_crits(binodals, uidx):

	raw_crits = np.array(binodals["groupings"][uidx]["raw_crits"])

	distances     = cdist(raw_crits, raw_crits)
	# Set distances between the same index to NaN
	np.fill_diagonal(distances, np.nan)

	# Get indices of upper triangle (excluding diagonal)
	upper_triangle_indices = np.triu_indices(raw_crits.shape[0], k=1)

	# Get distances corresponding to upper triangle indices
	distances = distances[upper_triangle_indices]

	print(distances)

	closeness = np.allclose(distances, distances[0], atol=1e-6)
	print(closeness)


	if closeness:
		return "cyclic"
	else:
		return "triumvirate"

#########################################
	
def line_check(crits, P):

	for i in range(len(crits)):
		for j in range(i+1, len(crits)):
			line = np.linspace(crits[i][0:2], crits[j][0:2], 500)
			to_check = ternary.stab_crit(line[:,0], line[:,1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) >= 0
			if not(to_check[10:-10].all()):
				return "triumvirate"
			

	return "triad"  # Concave hull computation failed, indicating non-convexity

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
					if np.sign(np.cross(central_axis[0:2], p1[0:2]-center[0:2])) == np.sign(np.cross(central_axis[0:2], p2[0:2]-center[0:2])):
						# print(f"Same side.")
						continue
					elif np.sign(np.cross(central_axis[0:2], p1[0:2]-center[0:2])) >=0 and np.sign(np.cross(central_axis[0:2], p2[0:2] - center[0:2])) <= 0:
						# print(f"Good!.")
						b_neg = np.vstack((b_neg, p2))
						b_pos = np.vstack((b_pos, p1))
					elif np.sign(np.cross(central_axis[0:2], p1[0:2]-center[0:2])) <=0 and np.sign(np.cross(central_axis[0:2], p2[0:2] - center[0:2])) >= 0:
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
			b_neg, b_pos     = solve_within(new_arm_neg, new_arm_pos, P, center, central_axis)
			arm_neg          = np.vstack((arm_neg, b_neg))
			arm_pos          = np.vstack((arm_pos, b_pos))
			arm_neg, arm_pos = clean_and_sort(arm_neg, arm_pos, center, central_axis)

		# find point with greatest gap
		arm_neg, keep = ternary.remove_close_rows(arm_neg, 1e-12)
		arm_pos       = arm_pos[keep]

		dist_neg      = np.linalg.norm(np.diff(arm_neg, axis=0), axis=1)
		max_dist_neg  = np.argmax(dist_neg)
		dist_pos      = np.linalg.norm(np.diff(arm_pos, axis=0), axis=1)
		max_dist_pos  = np.argmax(dist_pos)

		print(f"iterx = {iterx}", flush=True)
		print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
		iterx += 1
		if iterx > 100:
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

def populate_center(BINODALS, idx_tup, M=50):

	nv_c = BINODALS["crit_info"][BINODALS["groupings"][idx_tup]["center"  ]["idx"]]["norm_vec"]
	nv_p = BINODALS["crit_info"][BINODALS["groupings"][idx_tup]["positive"]["idx"]]["norm_vec"]
	nv_n = BINODALS["crit_info"][BINODALS["groupings"][idx_tup]["negative"]["idx"]]["norm_vec"]

	cp_c = P.crits[BINODALS["groupings"][idx_tup]["center"  ]["idx"]]
	cp_p = P.crits[BINODALS["groupings"][idx_tup]["positive"]["idx"]] 
	cp_n = P.crits[BINODALS["groupings"][idx_tup]["negative"]["idx"]] 

	Bneg     = BINODALS["groupings"][idx_tup]["negative"]["binodals"]
	Bpos     = BINODALS["groupings"][idx_tup]["positive"]["binodals"]
	Bc       = BINODALS["groupings"][idx_tup]["center"  ]["binodals"]

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
		BINODALS["groupings"][idx_tup]["negative"]["binodals"][0] = np.vstack((BINODALS["groupings"][idx_tup]["negative"]["binodals"][0], arm_neg))
		BINODALS["groupings"][idx_tup]["negative"]["binodals"][1] = np.vstack((BINODALS["groupings"][idx_tup]["negative"]["binodals"][1], arm_pos))
		BINODALS["groupings"][idx_tup]["negative"]["binodals"][0], BINODALS["groupings"][idx_tup]["negative"]["binodals"][1] = clean_and_sort(BINODALS["groupings"][idx_tup]["negative"]["binodals"][0], BINODALS["groupings"][idx_tup]["negative"]["binodals"][1], cp_n, nv_n) 

		d_neg = distance_from_axis(BINODALS["groupings"][idx_tup]["negative"]["binodals"][1][-1], cp_c, nv_c)
		
		# moving on to the positive side 
		to_insert_neg = np.linspace(Bpos[0][-1], Bneg[1][-1], M)
		to_insert_pos = np.linspace(Bpos[1][-1], Bc[1][0],    M)

		arm_neg, arm_pos = solve_central(to_insert_pos, to_insert_neg, cp_p, nv_p) 
		# only keep points to the negative of the central axis 
		adj_arm_neg = (cp_c[0:2]-arm_neg[:,0:2])/np.linalg.norm(cp_c[0:2]-arm_neg[:,0:2], axis=1).reshape(-1,1)
		signs       = np.sign(np.cross(nv_c, adj_arm_neg))
		arm_neg     = arm_neg[signs>=0]
		arm_pos     = arm_pos[signs>=0]
		BINODALS["groupings"][idx_tup]["positive"]["binodals"][0] = np.vstack((BINODALS["groupings"][idx_tup]["positive"]["binodals"][0], arm_neg))
		BINODALS["groupings"][idx_tup]["positive"]["binodals"][1] = np.vstack((BINODALS["groupings"][idx_tup]["positive"]["binodals"][1], arm_pos))
		BINODALS["groupings"][idx_tup]["positive"]["binodals"][0], BINODALS["groupings"][idx_tup]["positive"]["binodals"][1] = clean_and_sort(BINODALS["groupings"][idx_tup]["positive"]["binodals"][0], BINODALS["groupings"][idx_tup]["positive"]["binodals"][1], cp_p, nv_p) 

		d_pos = distance_from_axis(BINODALS["groupings"][idx_tup]["positive"]["binodals"][1][-1], cp_c, nv_c)
		iterx += 1
		if iterx > 10:
			print("Breaking out of the center populater...", flush=True)
			break
	
	return

##########################################

def solve_between_islands(arm_neg, arm_pos, hull_neg, hull_pos):

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
					if hull_neg.contains_point(p1[0:2]) and hull_pos.contains_point(p2[0:2]):
						# print(f"Good!")
						b_neg = np.vstack((b_neg, p1))
						b_pos = np.vstack((b_pos, p2))
					if hull_neg.contains_point(p2[0:2]) and hull_pos.contains_point(p1[0:2]):
						# print(f"Good!")
						b_neg = np.vstack((b_neg, p2))
						b_pos = np.vstack((b_pos, p1))
					else:
						# print(f"Something else...")
						continue
					break

	return b_neg, b_pos

##########################################

def populate_between_islands(BINODALS, idx_tup, M=50):

	print(f"idx_tup = {idx_tup}")

	arm_neg = BINODALS["island_scans"][idx_tup][idx_tup[0]]["binodal"]
	arm_pos = BINODALS["island_scans"][idx_tup][idx_tup[1]]["binodal"]

	print(f"arm_neg = {arm_neg[0:10]}",flush=True)
	print(f"arm_pos = {arm_pos[0:10]}",flush=True)

	# get the center 
	center = np.mean(np.vstack((BINODALS["island_scans"][idx_tup][idx_tup[0]]["binodal"], BINODALS["island_scans"][idx_tup][idx_tup[1]]["binodal"])), axis=0)

	dist_neg = np.linalg.norm(np.diff(arm_neg[:,0:2], axis=0), axis=1)
	dist_pos = np.linalg.norm(np.diff(arm_pos[:,0:2], axis=0), axis=1)

	max_dist_neg = np.argmax(dist_neg)
	max_dist_pos = np.argmax(dist_pos)
	
	iterx = 0
	while dist_neg[max_dist_neg] > 0.0005 or dist_pos[max_dist_pos] > 0.0005:
		
		if dist_neg[max_dist_neg] > dist_pos[max_dist_pos]:
			# treating the negative side 
			new_arm_neg = add_rows(arm_neg, M, max_dist_neg)
			new_arm_pos = add_rows(arm_pos, M, max_dist_neg)

			# solve the arms 
			b_neg, b_pos = solve_between_islands(new_arm_neg, new_arm_pos, BINODALS["stable_hulls"][idx_tup[0]], BINODALS["stable_hulls"][idx_tup[1]]) 

			# bring everything together
			arm_neg = np.vstack((arm_neg, b_neg))
			arm_pos = np.vstack((arm_pos, b_pos))

			# now sort everything
			angles  = np.arctan2(arm_neg[:,1]-center[1], arm_neg[:,0]-center[0])
			arm_neg = arm_neg[np.argsort(angles)]
			arm_pos = arm_pos[np.argsort(angles)]

			angles = angles[np.argsort(angles)]
			max_delta_theta = np.max(angles[1:] - angles[:-1])

			if max_delta_theta > np.pi/10:
				dists = np.linalg.norm(arm_neg[1:][0:2] - arm_neg[:-1][0:2], axis=1)
				max_dist_idx = np.argmax(dists)
				arm_neg = np.vstack((arm_neg[max_dist_idx+1:], arm_neg[:max_dist_idx+1]))
				arm_pos = np.vstack((arm_pos[max_dist_idx+1:], arm_pos[:max_dist_idx+1]))

		elif dist_neg[max_dist_neg] < dist_pos[max_dist_pos]:
			# treating the negative side 
			new_arm_neg = add_rows(arm_neg, M, max_dist_pos)
			new_arm_pos = add_rows(arm_pos, M, max_dist_pos)

			# solve the arms 
			b_neg, b_pos = solve_between_islands(new_arm_neg, new_arm_pos, BINODALS["stable_hulls"][idx_tup[0]], BINODALS["stable_hulls"][idx_tup[1]]) 

			# bring everything together
			arm_neg = np.vstack((arm_neg, b_neg))
			arm_pos = np.vstack((arm_pos, b_pos))

			# now sort everything
			angles  = np.arctan2(arm_pos[:,1]-center[1], arm_pos[:,0]-center[0])
			arm_neg = arm_neg[np.argsort(angles)]
			arm_pos = arm_pos[np.argsort(angles)]

			angles = angles[np.argsort(angles)]
			max_delta_theta = np.max(angles[1:] - angles[:-1])

			if max_delta_theta > np.pi/10:
				dists = np.linalg.norm(arm_neg[1:][0:2] - arm_neg[:-1][0:2], axis=1)
				max_dist_idx = np.argmax(dists)
				arm_neg = np.vstack((arm_neg[max_dist_idx+1:], arm_neg[:max_dist_idx+1]))
				arm_pos = np.vstack((arm_pos[max_dist_idx+1:], arm_pos[:max_dist_idx+1]))

		# find point with greatest gap
		arm_neg, keep = ternary.remove_close_rows(arm_neg, 1e-12)
		arm_pos       = arm_pos[keep]

		dist_neg      = np.linalg.norm(np.diff(arm_neg[:,0:2], axis=0), axis=1)
		max_dist_neg  = np.argmax(dist_neg)
		dist_pos      = np.linalg.norm(np.diff(arm_pos[:,0:2], axis=0), axis=1)
		max_dist_pos  = np.argmax(dist_pos)

		print(f"iterx = {iterx}", flush=True)
		print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
		iterx += 1
		if iterx > 20:
			print("Breaking out of the populater...", flush=True)
			break
	
	BINODALS["island_scans"][idx_tup][idx_tup[0]]["binodal"] = arm_neg
	BINODALS["island_scans"][idx_tup][idx_tup[1]]["binodal"] = arm_pos

	return

##########################################

def get_shared_idx(t1, t2):
	if t1[0] in t2:
		return t1[0]
	else:
		return t1[1]

##########################################

def triangle_finder(phi_, P):
	eq1 = P.sym_mu_ps.delta_mu_s(phi_[0], phi_[1], phi_[2], phi_[3])
	eq2 = P.sym_mu_ps.delta_mu_s(phi_[0], phi_[1], phi_[4], phi_[5])
	eq3 = P.sym_mu_ps.delta_mu_p(phi_[0], phi_[1], phi_[2], phi_[3])
	eq4 = P.sym_mu_ps.delta_mu_p(phi_[0], phi_[1], phi_[4], phi_[5])
	eq5 = P.sym_mu_ps.delta_mu_c(phi_[0], phi_[1], phi_[2], phi_[3])
	eq6 = P.sym_mu_ps.delta_mu_c(phi_[0], phi_[1], phi_[4], phi_[5])
	return [eq1, eq2, eq3, eq4, eq5, eq6]

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
	print(f"crits = {P.crits}", flush=True)
	print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	#=================================
	# create all the volume fractions
	p_s_space = np.arange(0.001, 1-0.001, 0.001)
	p_s = np.repeat(p_s_space, len(p_s_space))

	p_p = np.zeros(p_s.shape)
	for i in range (len(p_s_space)):
		p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))
	
	# curate all the volume fractions
	vals    = ternary.stab_crit(p_s, p_p, P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)
	to_keep = ~np.isnan(vals) & ~np.isinf(vals) & (vals > 0)

	p_s = p_s[to_keep]
	p_p = p_p[to_keep]
	phi = np.array([p_s, p_p]).T

	#=================================
	#=================================

	f = open(args.fb, 'rb')
	BINODALS = pickle.load(f)
	f.close()

	'''
	cols = ["black", "white", "grey", "pink", "limegreen", "gold"]
	for i in range(len(BINODALS["stable_hulls"])):
		for j in range(i+1, len(BINODALS["stable_hulls"])):
			keep_i = BINODALS["stable_hulls"][i].contains_points(BINODALS["island_scans"][(i,j)][i]["binodal"][:,0:2])
			keep_j = BINODALS["stable_hulls"][j].contains_points(BINODALS["island_scans"][(i,j)][j]["binodal"][:,0:2])
			keep   = np.logical_and(keep_i, keep_j)
			BINODALS["island_scans"][(i,j)][i]["binodal"] =  BINODALS["island_scans"][(i,j)][i]["binodal"][keep]
			BINODALS["island_scans"][(i,j)][j]["binodal"] =  BINODALS["island_scans"][(i,j)][j]["binodal"][keep]
	'''

	for ikey, key in enumerate(BINODALS["hull_info"]["binodal"]):
		# print(key)
		if key[-1]=="two_phase":
			B1 = key[0][0]
			B2 = key[0][1]
			if len(B1) == 0:
				continue
			B1, keep = ternary.remove_close_rows(B1, 1e-6)
			B2       = B2[keep]

			# c1  = P.crits[0]
			# c2  = P.crits[1]
			# mid = (c1+c2)/2
			# mid_axis    = (mid - c1)[0:2]/np.linalg.norm((mid-c1)[0:2])
			# adj_neg_arm = (B1[:,0:2]-mid[0:2])/np.linalg.norm(B1[:,0:2]-mid[0:2], axis=1)[:,np.newaxis]
			# angles      = np.arccos(np.sum(adj_neg_arm * mid_axis, axis=1))
			# B1  = B1[np.argsort(angles)]
			# B2  = B2[np.argsort(angles)]

			ax.plot(B1[:,0], 1-B1[:,0]-B1[:,1], B1[:,1], c="white", lw=1, label="binodal", zorder=10)
			ax.plot(B2[:,0], 1-B2[:,0]-B2[:,1], B2[:,1], c="black", lw=1, label="binodal", zorder=10)

			# ax.scatter(B1[:,0], 1-B1[:,0]-B1[:,1], B1[:,1], c="white", s=1, label="binodal", zorder=10)
			# ax.scatter(B2[:,0], 1-B2[:,0]-B2[:,1], B2[:,1], c="black", s=1, label="binodal", zorder=10)

			for i in range(0, len(B1)):
				ax.plot([B1[i,0],B2[i,0]], [1-B1[i,0]-B1[i,1],1-B2[i,0]-B2[i,1]], [B1[i,1], B2[i,1]], lw=1, ls='--', c="pink", zorder=5)

		else:
			B = key[0]
			ax.plot(np.hstack((B[:,0],B[0,0])), np.hstack((1-B[:,0]-B[:,1],1-B[0,0]-B[0,1])), np.hstack((B[:,1],B[0,1])), lw=1, c="black", zorder=10)
	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
	
	exit()
	

	'''
	# process the intersections
	keys = list(BINODALS["island_scans"].keys())
	for i in range(len(keys)):
		if not ("intersections" in list(BINODALS["island_scans"][keys[i]].keys())):
			BINODALS["island_scans"][keys[i]]["intersections"] = dict()
		for j in range(i+1, len(keys)):
			if not ("intersections" in list(BINODALS["island_scans"][keys[j]].keys())):
				BINODALS["island_scans"][keys[j]]["intersections"] = dict()
			if keys[i][0] in keys[j]:
				BINODALS["island_scans"][keys[i]]["intersections"][keys[j]] = dict()
				BINODALS["island_scans"][keys[j]]["intersections"][keys[i]] = dict() 

				idx_in_j = keys[j].index(keys[i][0])
				L1 = LineString(BINODALS["island_scans"][keys[i]][keys[i][0]       ]["binodal"][:, 0:2])
				L2 = LineString(BINODALS["island_scans"][keys[j]][keys[j][idx_in_j]]["binodal"][:, 0:2])

				intersection = L1.intersection(L2)
				if intersection.is_empty:
					BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = False
					BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = False
					print(f"No intersection.", flush=True)
				else:
					try:
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = True
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = True
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["poi"]        = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["poi"]        = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
					except AttributeError: 
						print(f"Error in keys[i] = {keys[i]}, and keys[j] = {keys[j]}", flush=True)
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = False
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = False
					# BINODALS["island_scans"][keys[i]][keys[i][0]]["intersections"][keys[j]] = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
					# BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]             = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
			elif keys[i][1] in keys[j]:
				BINODALS["island_scans"][keys[i]]["intersections"][keys[j]] = dict()
				BINODALS["island_scans"][keys[j]]["intersections"][keys[i]] = dict()

				idx_in_j = keys[j].index(keys[i][1])
				L1 = LineString(BINODALS["island_scans"][keys[i]][keys[i][1]       ]["binodal"][:, 0:2])
				L2 = LineString(BINODALS["island_scans"][keys[j]][keys[j][idx_in_j]]["binodal"][:, 0:2])

				intersection = L1.intersection(L2)
				if intersection.is_empty:
					BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = False
					BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = False
				else:
					try:
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = True
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = True
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["poi"]        = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["poi"]        = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
					except AttributeError:
						print(f"Error in keys[i] = {keys[i]}, and keys[j] = {keys[j]}", flush=True)
						BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = False
						BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = False
			else:
				# no common island, no intersections 
				print(f"No common island, no intersections. Jump past...", flush=True)
				pass
	
	# now, let's clip the binodals.
	keys = list(BINODALS["island_scans"])
	triangle_set  = set()

	# get the sets of points that make the triangles
	for key in BINODALS["island_scans"]:
		kcount = 0
		triangle = [key]
		for k2 in BINODALS["island_scans"][key]["intersections"]:
			print(f"key = {key}, k2 = {k2}...", end=' ', flush=True)
			print(f'First bool = {BINODALS["island_scans"][key]["intersections"][k2]["boolean"]}, second bool = {BINODALS["island_scans"][k2]["intersections"][key]["boolean"]}...', end=' ', flush=True)
			if BINODALS["island_scans"][key]["intersections"][k2]["boolean"] and BINODALS["island_scans"][k2]["intersections"][key]["boolean"]:
				print(f"Appending.")
				triangle.append(k2)
			else:
				print(f"Not appending.")

		if len(triangle) == 3:
			print(f"triangle = {triangle}")
			print(f"triangle[0] = {triangle[0]}, triangle[1] = {triangle[1]}, triangle[2] = {triangle[2]}")
			if BINODALS["island_scans"][triangle[0]]["intersections"][triangle[1]]["boolean"] \
				and BINODALS["island_scans"][triangle[0]]["intersections"][triangle[2]]["boolean"] \
				and BINODALS["island_scans"][triangle[1]]["intersections"][triangle[2]]["boolean"]:
				triangle.sort()
				triangle = tuple(triangle)
				triangle_set.add(triangle)
			else:
				pass


	print(triangle_set)

	# go to each triangle, and start finding the three points
	triangular_hulls = []
	for triangle_keys in triangle_set:
		# get the three points 
		triangle_points = list() 
		for i in range(len(triangle_keys)):
			for j in range(i+1, len(triangle_keys)):
				triangle_points.append(BINODALS["island_scans"][triangle_keys[i]]["intersections"][triangle_keys[j]]["poi"][0:2])
		
		root = fsolve(triangle_finder, [triangle_points[0][0], triangle_points[0][1], triangle_points[1][0], triangle_points[1][1], triangle_points[2][0], triangle_points[2][1]], args=(P))

		print(f"root = {root}", flush=True)
		closeness = np.isclose(triangle_finder(root, P), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
		print(f"closeness = {closeness}", flush=True)
		if closeness.all():
			print("Worked!")
		else: 
			print("This ain't a triangle...", flush=True)

		triangular_hull = np.array([[root[0], root[1]], [root[2], root[3]], [root[4], root[5]]])

		hull      = ConvexHull(triangular_hull[:,0:2])
		hull_path = Path(triangular_hull[:,0:2][hull.vertices])
		triangular_hulls.append([triangular_hull, hull_path])
		# BINODALS["hull_info"]["triangles"].append(hull_path)
		# BINODALS["hull_info"]["function"].append(hull_path)
		# BINODALS["hull_info"]["binodal"].append([triangular_hull, "three_phase"])

	to_pop = []
	for i in range(len(triangular_hulls)):
		for j in range(len(triangular_hulls)):
			if i == j:
				continue 
			if (triangular_hulls[i][1].contains_points(triangular_hulls[j][0])).any():
				to_pop.append(j)

	for i in range(len(triangular_hulls)):
		if not(i in to_pop):
			BINODALS["hull_info"]["triangles"].append(triangular_hulls[i][1])
			BINODALS["hull_info"]["function"].append(triangular_hulls[i][1])
			BINODALS["hull_info"]["binodal"].append([triangular_hulls[i][0], "three_phase"])

	# go to each BINODAL and plug them into hull_info
	for key in BINODALS["island_scans"]:
		print(f"key = {key}", flush=True)
		BINODALS["hull_info"]["binodal"].append( [[BINODALS["island_scans"][key][key[0]]["binodal"], BINODALS["island_scans"][key][key[1]]["binodal"]], "two_phase"] )
		if BINODALS["island_scans"][key][key[0]]["binodal"].shape[0] == 0 and BINODALS["island_scans"][key][key[1]]["binodal"].shape[0] == 0:
			BINODALS["hull_info"]["function"].append(None)
		else:	
			combined = np.vstack(([BINODALS["island_scans"][key][key[0]]["binodal"], BINODALS["island_scans"][key][key[1]]["binodal"]]))
			hull = ConvexHull(combined[:,0:2])
			hull_path = Path(combined[:,0:2][hull.vertices])
			BINODALS["hull_info"]["function"].append(hull_path)

	for idx, H in enumerate(BINODALS["hull_info"]["binodal"]):
		if H[-1] == "two_phase" and H[0][0].shape[0] !=0 and H[0][1].shape[1] != 0:
			to_keep = list()
			for i in range(len(H[0][0])):
				check = True
				line = np.linspace(H[0][0][i][0:2], H[0][1][i][0:2], 1000)
				for hull_path in BINODALS["hull_info"]["triangles"]:
					if (hull_path.contains_points(line)).any():
						check = False
						break
					else:
						continue
				if check:
					to_keep.append(i)
			BINODALS["hull_info"]["binodal"][idx][0][0] = BINODALS["hull_info"]["binodal"][idx][0][0][to_keep]
			BINODALS["hull_info"]["binodal"][idx][0][1] = BINODALS["hull_info"]["binodal"][idx][0][1][to_keep]

			if BINODALS["hull_info"]["binodal"][idx][0][0].shape[0] < 3 or BINODALS["hull_info"]["binodal"][idx][0][1].shape[0] < 3:
				BINODALS["hull_info"]["function"][idx] = None
			else:
				combined  = np.vstack((BINODALS["hull_info"]["binodal"][idx][0][0], BINODALS["hull_info"]["binodal"][idx][0][1]))
				hull      = ConvexHull(combined[:,0:2])
				hull_path = Path(combined[:,0:2][hull.vertices])
				BINODALS["hull_info"]["function"][idx] = hull_path

		else:
			continue

	for idx, H in enumerate(BINODALS["hull_info"]["binodal"]):
		if H[0][0].shape[0] == 0 and H[0][1].shape[0] == 0:
			continue
		elif H[-1] == "two_phase":
			ax.scatter(H[0][0][:,0], 1-H[0][0][:,0]-H[0][0][:,1], H[0][0][:,1], s=0.5, c='black')
			ax.scatter(H[0][1][:,0], 1-H[0][1][:,0]-H[0][1][:,1], H[0][1][:,1], s=0.5, c='white')
		elif H[-1] == "three_phase":
			ax.plot(np.hstack([H[0][:,0],H[0][0,0]]),\
			np.hstack([1-H[0][:,0]-H[0][:,1], 1-H[0][0,0]-H[0][0,1]]), np.hstack([H[0][:,1], H[0][0,1]]), c='slategray', lw=1)


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
		fig.savefig (f"probe-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	
	print(f"done!", flush=True)
	'''
	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

