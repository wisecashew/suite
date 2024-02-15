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
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

import argparse
parser = argparse.ArgumentParser(description='This will calculate a binodal.')
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
parser.add_argument('--mesh-pkl',            metavar='MPKL',    dest='mpkl',          type=str,   action='store', help='extract information about the mesh (default: None).',                                   default=None)
parser.add_argument('--crit-pkl',            metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                                default=None)
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

def transform_islands(islands, mesh):
	hull_paths = []
	for idx, island in enumerate(islands):
		phi_s = mesh[0][islands[idx][:,0], islands[idx][:,1]]
		phi_p = mesh[1][islands[idx][:,0], islands[idx][:,1]]
		islands[idx] = np.array([phi_s, phi_p]).T
		# islands[idx] = np.array(np.vstack([0.001+(0.999-0.001)/args.sd*island[:,1], 0.001+(0.999-0.001)/args.sd*island[:,0]])).T 
		hull = ConvexHull(islands[idx])
		hull_paths.append(Path(islands[idx][hull.vertices]))
		# print(hull_paths[-1])
		# print(hull_paths[-1].contains_point(p1[0:2]))
		# print(hull_paths[-1].contains_point(p2[0:2]))

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

	while dist_neg[max_dist_neg] > 0.001 or dist_pos[max_dist_pos] > 0.001:

		# print(f"Negative points are: {arm_neg[max_dist_neg], arm_neg[max_dist_neg+1]}")
		# print(f"Positive points are: {arm_pos[max_dist_pos], arm_pos[max_dist_pos+1]}")

		if dist_neg[max_dist_neg] > dist_pos[max_dist_pos] and dist_neg[max_dist_neg] > 0.001:
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


		elif dist_pos[max_dist_pos] > dist_neg[max_dist_neg] and dist_neg[max_dist_neg] > 0.001:
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
		arm_neg, keep = ternary.remove_close_rows(arm_neg, 1e-15)
		arm_pos       = arm_pos[keep]

		dist_neg      = np.linalg.norm(np.diff(arm_neg, axis=0), axis=1)
		max_dist_neg  = np.argmax(dist_neg)
		dist_pos      = np.linalg.norm(np.diff(arm_pos, axis=0), axis=1)
		max_dist_pos  = np.argmax(dist_pos)

		print(f"iterx = {iterx}", flush=True)
		print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
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
					# print(hull_neg)
					# print(hull_pos)
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
				dists = np.linalg.norm(arm_pos[1:][0:2] - arm_pos[:-1][0:2], axis=1)
				max_dist_idx = np.argmax(dists)
				arm_neg = np.vstack((arm_neg[max_dist_idx+1:], arm_neg[:max_dist_idx+1]))
				arm_pos = np.vstack((arm_pos[max_dist_idx+1:], arm_pos[:max_dist_idx+1]))

		# find point with greatest gap
		arm_neg, keep = ternary.remove_close_rows(arm_neg, 1e-15)
		arm_pos       = arm_pos[keep]

		dist_neg      = np.linalg.norm(np.diff(arm_neg[:,0:2], axis=0), axis=1)
		max_dist_neg  = np.argmax(dist_neg)
		dist_pos      = np.linalg.norm(np.diff(arm_pos[:,0:2], axis=0), axis=1)
		max_dist_pos  = np.argmax(dist_pos)

		print(f"iterx = {iterx}", flush=True)
		print(f"maximum gap in dist_neg = {dist_neg[max_dist_neg]}, maximum gap in dist_pos = {dist_pos[max_dist_pos]}", flush=True)
		iterx += 1
		if iterx > 5: # 50:
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

def extend_binodals(BINODALS, P, idx_tup):

	# Determine extension direction based on the difference between the last two points
	neg_arm = BINODALS["island_scans"][idx_tup][idx_tup[0]]["binodal"]
	pos_arm = BINODALS["island_scans"][idx_tup][idx_tup[1]]["binodal"]

	# get an extension direction 
	neg_dir = (neg_arm[-1][0:2]-neg_arm[-2][0:2])/np.linalg.norm(neg_arm[-1][0:2]-neg_arm[-2][0:2])
	pos_dir = (pos_arm[-1][0:2]-pos_arm[-2][0:2])/np.linalg.norm(pos_arm[-1][0:2]-pos_arm[-2][0:2])

	# define iterations for the while loops
	forw_extension_maxiter = 1e+4 
	back_extension_maxiter = 1e+4

	# define certain parameters for the loop
	forw_extension_iter = 0
	forw_scaler = 0.05
	while forw_extension_iter < forw_extension_maxiter:
		forw_extension_iter += 1

		neg_ext = neg_arm[-1][0:2] + forw_scaler*neg_dir 
		pos_ext = pos_arm[-1][0:2] + forw_scaler*pos_dir 

		def mu_equations(phi):
			eq1 = P.sym_mu_ps.delta_mu_s(neg_ext[0], phi[0], phi[1], phi[2]) 
			eq2 = P.sym_mu_ps.delta_mu_p(neg_ext[0], phi[0], phi[1], phi[2]) 
			eq3 = P.sym_mu_ps.delta_mu_c(neg_ext[0], phi[0], phi[1], phi[2]) 
			return [eq1, eq2, eq3]	

		root = fsolve(mu_equations, [neg_ext[1], pos_ext[0], pos_ext[1]])
		if (np.abs(np.array(mu_equations(root)))>1e-6).any():
			forw_scaler = forw_scaler/2
		elif ternary.stab_crit(neg_ext[0], neg_ext[1], P.vs, P.vc, P.vp, P.chi_ps,  P.chi_pc, P.chi_sc) < 0 or ternary.stab_crit(pos_ext[0], pos_ext[1], P.vs, P.vc, P.vp, P.chi_ps,  P.chi_pc, P.chi_sc) < 0:
			forw_scaler = forw_scaler/2
		else:
			# print(f"Found a valid extension on the front end!")
			p_neg   = np.array([neg_ext[0], root[0], 1-neg_ext[0]-root[0]])
			p_pos   = np.array([root[1], root[2], 1-root[1]-root[2]])
			neg_arm = np.vstack((neg_arm, p_neg))
			pos_arm = np.vstack((pos_arm, p_pos))
			continue

	#################
	
	neg_dir = (neg_arm[0][0:2]-neg_arm[1][0:2])/np.linalg.norm(neg_arm[0][0:2]-neg_arm[1][0:2])
	pos_dir = (pos_arm[0][0:2]-pos_arm[1][0:2])/np.linalg.norm(pos_arm[0][0:2]-pos_arm[1][0:2])

	back_extension_iter = 0
	back_scaler         = 0.05
	while back_extension_iter < back_extension_maxiter:
		back_extension_iter += 1
		
		neg_ext = neg_arm[0][0:2] + back_scaler*neg_dir
		pos_ext = pos_arm[0][0:2] + back_scaler*pos_dir

		def mu_equations(phi):
			eq1 = P.sym_mu_ps.delta_mu_s(neg_ext[0], phi[0], phi[1], phi[2]) 
			eq2 = P.sym_mu_ps.delta_mu_p(neg_ext[0], phi[0], phi[1], phi[2]) 
			eq3 = P.sym_mu_ps.delta_mu_c(neg_ext[0], phi[0], phi[1], phi[2]) 
			return [eq1, eq2, eq3]	

		root = fsolve(mu_equations, [neg_ext[1], pos_ext[0], pos_ext[1]])
		if (np.abs(np.array(mu_equations(root)))>1e-6).any():
			back_scaler = back_scaler/2
		elif ternary.stab_crit(neg_ext[0], neg_ext[1], P.vs, P.vc, P.vp, P.chi_ps,  P.chi_pc, P.chi_sc) < 0 or ternary.stab_crit(pos_ext[0], pos_ext[1], P.vs, P.vc, P.vp, P.chi_ps,  P.chi_pc, P.chi_sc) < 0:
			back_scaler = back_scaler/2
		else:
			# print(f"Found a valid extension on the back end!")
			p_neg   = np.array([neg_ext[0], root[0], 1-neg_ext[0]-root[0]])
			p_pos   = np.array([root[1], root[2], 1-root[1]-root[2]])
			neg_arm = np.vstack((p_neg, neg_arm))
			pos_arm = np.vstack((p_pos, pos_arm))


	BINODALS["island_scans"][idx_tup][idx_tup[0]]["binodal"] = neg_arm 
	BINODALS["island_scans"][idx_tup][idx_tup[1]]["binodal"] = pos_arm

	return 


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
	# fig.savefig("to_del.png", dpi=1200)
	# exit()
	# P.tangent_tracing(ax)
	# print("Plotted out the tangent trace!", flush=True)
	
	# extract island information
	f = open(args.spkl, 'rb')
	stable_islands = pickle.load(f)
	f.close()

	f = open(args.upkl, 'rb')
	unstable_islands = pickle.load(f)
	f.close()

	f = open(args.mpkl, 'rb')
	mesh = pickle.load(f)
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
	#=================================
	#=================================
	# calculate the stable centers in your system
	stable_centers   = []
	hull_paths_s = transform_islands(stable_islands, mesh)

	cols = ["limegreen", "darkred", "steelblue", "lavender"]
	for sidx, si in enumerate(stable_islands):
		stable_centers.append(np.mean(si, axis=0))
		# ax.scatter(si[:,0],1-si[:,0]-si[:,1], si[:,1], c=cols[sidx])
		# print(f"phi_s = {np.min(si[:,0])}, phi_p = {np.min(si[:,1])}, phi_c = {np.min(1-si[:,0]-si[:,1])}")

	unstable_centers = []
	hull_paths_u = transform_islands(unstable_islands, mesh)
	
	cols = ["white", "black", "slategray"]
	for uidx, ui in enumerate(unstable_islands):
		unstable_centers.append(np.mean(ui, axis=0))
		# ax.scatter(ui[:,0],1-ui[:,0]-ui[:,1], ui[:,1], c=cols[uidx], s=0.5)
		# print(f"col = {cols[uidx]}")
		# print(f"phi_s = {np.min(ui[:,0])}, phi_p = {np.min(ui[:,1])}, phi_c = {np.min(1-ui[:,0]-ui[:,1])}")
	unstable_centers = np.array(unstable_centers)
	
	# fig.savefig("testing_mesh.png", dpi=1200, bbox_inches="tight")
	# exit()
	# ax.scatter(stable_centers[:,0],   1-stable_centers[:,0]-stable_centers[:,1],       stable_centers[:,1], c='hotpink', s=2, zorder=200)
	# ax.scatter(unstable_centers[:,0], 1-unstable_centers[:,0]-unstable_centers[:,1], unstable_centers[:,1], c='orange' , s=2, zorder=200)
	
	#=================================
	def triangle_finder(phi_):
		eq1 = P.sym_mu_ps.delta_mu_s(phi_[0], phi_[1], phi_[2], phi_[3])
		eq2 = P.sym_mu_ps.delta_mu_s(phi_[0], phi_[1], phi_[4], phi_[5])
		eq3 = P.sym_mu_ps.delta_mu_p(phi_[0], phi_[1], phi_[2], phi_[3])
		eq4 = P.sym_mu_ps.delta_mu_p(phi_[0], phi_[1], phi_[4], phi_[5])
		eq5 = P.sym_mu_ps.delta_mu_c(phi_[0], phi_[1], phi_[2], phi_[3])
		eq6 = P.sym_mu_ps.delta_mu_c(phi_[0], phi_[1], phi_[4], phi_[5])
		return [eq1, eq2, eq3, eq4, eq5, eq6]
	#=================================
	
	BINODALS = dict()
	BINODALS["stable_hulls"]     = hull_paths_s 
	BINODALS["unstable_hulls"]   = hull_paths_u
	BINODALS["stable_islands"]   = stable_islands
	BINODALS["unstable_islands"] = unstable_islands
	BINODALS["groupings"]      = dict()
	BINODALS["island_scans"]   = dict()
	BINODALS["crit_info"]      = dict()
	BINODALS["hull_info"]      = dict()
	BINODALS["hull_info"]["triangles"] = list()
	BINODALS["hull_info"]["function"]  = list()
	BINODALS["hull_info"]["binodal" ]  = list()  

	# compile together all the geometric information of the critical points
	# theta = np.linspace(0, 2*np.pi, 1000)
	for cidx, crit in enumerate(P.crits):
		BINODALS["crit_info"][cidx] = dict()
		BINODALS["crit_info"][cidx]["tang_slope"] = tangent.tangent2(P.vs, P.vc, P.vp, crit[0], crit[1], P.chi_pc, P.chi_ps, P.chi_sc, P.spinodal.root_up_s, P.spinodal.root_lo_s)
		BINODALS["crit_info"][cidx]["tang_vec"]   = np.array([1, BINODALS["crit_info"][cidx]["tang_slope"]])/np.sqrt(1+BINODALS["crit_info"][cidx]["tang_slope"]**2)
		BINODALS["crit_info"][cidx]["norm_slope"] = -1/BINODALS["crit_info"][cidx]["tang_slope"]
		norm_vec = np.array([1, BINODALS["crit_info"][cidx]["norm_slope"]])/np.sqrt(1+BINODALS["crit_info"][cidx]["norm_slope"]**2)
		test_point = crit[0:2] + 0.01*norm_vec 
		if ternary.stab_crit(test_point[0], test_point[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) > 0:
			BINODALS["crit_info"][cidx]["norm_vec"]   = norm_vec
		else:
			BINODALS["crit_info"][cidx]["norm_vec"]   = -norm_vec

		print(f"crit = {crit}")
		'''
		L = np.linspace(0,0.01,500) 
		N = crit[0:2] - BINODALS["crit_info"][cidx]["norm_vec"][0:2]*L[:, np.newaxis]
		# print(f"N = {N}")
		ax.scatter(N[:,0], 1-N[:,0]-N[:,1], N[:,1], c='pink', s=0.5)
		for idx, hpu in enumerate(hull_paths_u):
			# ax.scatter(unstable_islands[idx][:,0], 1-unstable_islands[idx][:,0]-unstable_islands[idx][:,1], unstable_islands[idx][:,1], c='white', s=0.5)
			if (hpu.contains_points(N)).any():
				uidx = idx 
				break 
		
		N = crit[0:2] + BINODALS["crit_info"][cidx]["norm_vec"][0:2] *L[:, np.newaxis]
		ax.scatter(N[:,0], 1-N[:,0]-N[:,1], N[:,1], c='darkred', s=0.5)
		# print(f"N = {N}")
		for idx, hps in enumerate(hull_paths_s):
			# ax.scatter(stable_islands[idx][:,0], 1-stable_islands[idx][:,0]-stable_islands[idx][:,1], stable_islands[idx][:,1], c='slategray', s=0.5)
			if (hps.contains_points(N)).any():
				sidx = idx
				break
		'''
		min_dists_unstable = []
		for idx, ui in enumerate(BINODALS["unstable_islands"]):
			min_dist = np.min(np.linalg.norm(ui[:,0:2] - crit[0:2], axis=1))
			min_dists_unstable.append(min_dist)
		uidx = np.argmin(min_dists_unstable)

		min_dists_stable = []
		for idx, si in enumerate(BINODALS["stable_islands"]):
			min_dist = np.min(np.linalg.norm(si[:,0:2] - crit[0:2], axis=1))
			min_dists_stable.append(min_dist)
		sidx = np.argmin(min_dists_stable)

		print(f"uidx, sidx = {uidx, sidx}")
		
		# uidx = np.argmin(np.linalg.norm(unstable_centers - crit[0:2], axis=1))
		# sidx = np.argmin(np.linalg.norm(stable_centers   - crit[0:2], axis=1))

		if (uidx,sidx) in list(BINODALS["groupings"].keys()):
			BINODALS["groupings"][(uidx, sidx)]["raw_list" ].append(cidx)
			BINODALS["groupings"][(uidx, sidx)]["raw_crits"].append(crit)
		else:
			BINODALS["groupings"][(uidx, sidx)] = dict()
			BINODALS["groupings"][(uidx, sidx)]["binodals"]  = []
			BINODALS["groupings"][(uidx, sidx)]["raw_list"]  = [cidx] 
			BINODALS["groupings"][(uidx, sidx)]["raw_crits"] = [crit]

	# fig.savefig("to_del.png", dpi=1200, bbox_inches="tight")
	# exit()

	# time to order the groupings
	for idx_tup in list(BINODALS["groupings"].keys()):
		
		if len(BINODALS["groupings"][idx_tup]["raw_list"]) == 1:
			BINODALS["groupings"][idx_tup]["identity"]     = "unity"
			BINODALS["groupings"][idx_tup]["alpha"]        = dict()
			BINODALS["groupings"][idx_tup]["alpha"]["idx"] = BINODALS["groupings"][idx_tup]["raw_list"][0]

		elif len(BINODALS["groupings"][idx_tup]["raw_list"]) == 2:
			BINODALS["groupings"][idx_tup]["identity"]     = "dyad"
			BINODALS["groupings"][idx_tup]["alpha"]        = dict()
			BINODALS["groupings"][idx_tup]["beta" ]        = dict()
			BINODALS["groupings"][idx_tup]["alpha"]["idx"] = BINODALS["groupings"][idx_tup]["raw_list"][0]
			BINODALS["groupings"][idx_tup]["beta" ]["idx"] = BINODALS["groupings"][idx_tup]["raw_list"][1]

		elif len(BINODALS["groupings"][idx_tup]["raw_list"]) == 3:
			print(f"Performing a line check...", end= ' ', flush=True)
			BINODALS["groupings"][idx_tup]["identity"] = line_check(BINODALS["groupings"][idx_tup]["raw_crits"], P)
			print("done!", flush=True)

			if BINODALS["groupings"][idx_tup]["identity"] == "triad":
				titles = ["alpha", "beta", "gamma"]
				for ridx, raw_idx in enumerate(BINODALS["groupings"][idx_tup]["raw_list"]):
					BINODALS["groupings"][idx_tup][titles[ridx]] = dict()
					BINODALS["groupings"][idx_tup][titles[ridx]]["idx"] = raw_idx 

			elif BINODALS["groupings"][idx_tup]["identity"] == "triumvirate":
				BINODALS["groupings"][idx_tup]["center"  ] = dict()
				BINODALS["groupings"][idx_tup]["positive"] = dict()
				BINODALS["groupings"][idx_tup]["negative"] = dict()

				# find the middle point
				com = np.mean(BINODALS["groupings"][idx_tup]["raw_crits"], axis=0)
				BINODALS["groupings"][idx_tup]["center"]["idx"] = BINODALS["groupings"][idx_tup]["raw_list"][np.argmin(np.linalg.norm(np.array(BINODALS["groupings"][idx_tup]["raw_crits"])[:,0:2]- com[0:2], axis=1))]
				central_axis = BINODALS["crit_info"][BINODALS["groupings"][idx_tup]["center"]["idx"]]["norm_vec"]

				for c in BINODALS["groupings"][idx_tup]["raw_list"]:
					if c == BINODALS["groupings"][idx_tup]["center"]["idx"]:
						continue 
					else:
						deviation = (P.crits[c][0:2] - P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2])/np.linalg.norm(P.crits[c][0:2] - P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2])
						sign = np.sign(np.cross(BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["center"]["idx"] ]["norm_vec"][0:2], deviation[0:2]))
						if sign > 0:
							BINODALS["groupings"][idx_tup]["positive"]["idx"] = c
						else:
							BINODALS["groupings"][idx_tup]["negative"]["idx"] = c

	# end of setup. Now to move into specifics. 
	# now go through each idx_tup to start getting the binodals using normal pushing
	for idx_tup in list(BINODALS["groupings"].keys()):
		print(f"idx_tup = {idx_tup}", flush=True)

		if BINODALS["groupings"][idx_tup]["identity"] == "unity":
			P.tangent_tracing_unity(BINODALS, idx_tup)
			norm_vec = BINODALS["crit_info"][BINODALS["groupings"][idx_tup]["alpha"]["idx"]]["norm_vec"]
			center   = P.crits[BINODALS["groupings"][idx_tup]["alpha"]["idx"]]
			neg_arm  = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0]
			pos_arm  = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1]
			neg_sol, pos_sol = clean_and_sort(neg_arm, pos_arm, center, norm_vec)
			neg_sol = np.vstack((center[0:2], neg_sol))
			pos_sol = np.vstack((center[0:2], pos_sol))
			mask_neg = (ternary.stab_crit(neg_sol[:,0], neg_sol[:,1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)>=0)
			mask_pos = (ternary.stab_crit(pos_sol[:,0], pos_sol[:,1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)>=0)
			neg_sol = neg_sol[mask_neg]
			pos_sol = pos_sol[mask_pos]
			print(f"neg_sol = {neg_sol[-10:,0], neg_sol[-10:,1], 1-neg_sol[-10:,0]-neg_sol[-10:,1]}")
			print(f"pos_sol = {pos_sol[-10:,0], pos_sol[-10:,1], 1-pos_sol[-10:,0]-pos_sol[-10:,1]}")
			if len(neg_sol) == 0 or len(pos_sol) == 0:
				BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0] = np.empty((0,3))
				BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1] = np.empty((0,3))
				BINODALS["hull_info"]["function"].append(None)
				BINODALS["hull_info"]["binodal"].append([[np.empty((0,3)), np.empty((0,3))]])	
			else:
				combined  = np.vstack((neg_sol, pos_sol))
				hull      = ConvexHull(combined[:,0:2])
				hull_path = Path(combined[:,0:2][hull.vertices])
				BINODALS["hull_info"]["function"].append(hull_path) 
				BINODALS["hull_info"]["binodal"].append([BINODALS["groupings"][idx_tup]["alpha"]["binodals"], "two_phase"])

		elif BINODALS["groupings"][idx_tup]["identity"] == "dyad":
			print(f"Entering the tangent trace...", flush=True)
			P.tangent_tracing_dyad(BINODALS, idx_tup)
			print(f"done with tangent tracing!", flush=True)

			other_center = P.crits[BINODALS["groupings"][idx_tup]["raw_list"][0]]
			center = P.crits[BINODALS["groupings"][idx_tup]["raw_list"][1]]
			axis   = (P.crits[BINODALS["groupings"][idx_tup]["raw_list"][0]] - P.crits[BINODALS["groupings"][idx_tup]["raw_list"][1]])/np.linalg.norm(P.crits[BINODALS["groupings"][idx_tup]["raw_list"][0]]-P.crits[BINODALS["groupings"][idx_tup]["raw_list"][1]])
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], keep = ternary.remove_close_rows(BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], 1e-15)
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1]       = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1][keep]
			middle   = (center+other_center)/2
			mid_axis = (middle - center)[0:2]/np.linalg.norm((middle-center)[0:2])
			adj_neg_arm = (BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][:,0:2]-middle[0:2])/np.linalg.norm(BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][:,0:2]-middle[0:2]).reshape(-1,1)
			angles = np.arccos( np.sum(adj_neg_arm * mid_axis, axis=1) )
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0] = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][np.argsort(angles)]
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1] = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1][np.argsort(angles)]

			print(f"Entering add_and_solve...", flush=True)
			neg_arm, pos_arm = add_and_solve(BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1], P, center, mid_axis, 50)
			print(f"Got out of add and solve!", flush=True)

			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], neg_arm))
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1], pos_arm))
			adj_neg_arm = (BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][:,0:2]-middle[0:2])/np.linalg.norm(BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][:,0:2]-middle[0:2]).reshape(-1,1)
			angles = np.arccos( np.sum(adj_neg_arm * mid_axis, axis=1) )
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0] = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0][np.argsort(angles)]
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1] = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1][np.argsort(angles)]
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], keep = ternary.remove_close_rows(BINODALS["groupings"][idx_tup]["alpha"]["binodals"][0], 1e-15)
			BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1]       = BINODALS["groupings"][idx_tup]["alpha"]["binodals"][1][keep]

			combined  = np.vstack(BINODALS["groupings"][idx_tup]["alpha"]["binodals"])
			if (ternary.stab_crit(combined[:,0], combined[:,1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)<-1e-6).any():
				BINODALS["hull_info"]["function"].append(None) 
				BINODALS["hull_info"]["binodal"].append([np.empty((0,3)), np.empty((0,3))])
			else:
				hull      = ConvexHull(combined[:,0:2])
				hull_path = Path(combined[:,0:2][hull.vertices])
				BINODALS["hull_info"]["function"].append(hull_path) 
				BINODALS["hull_info"]["binodal"].append([BINODALS["groupings"][idx_tup]["alpha"]["binodals"], "two_phase"])
			# B = BINODALS["groupings"][idx_tup]["alpha"]["binodals"]
			# if args.pb:
			# 	ax.scatter(B[0][:,0], 1-B[0][:,0]-B[0][:,1], B[0][:,1], s=0.5, c='black')
			# 	ax.scatter(B[1][:,0], 1-B[1][:,0]-B[1][:,1], B[1][:,1], s=0.5, c='white')

		elif BINODALS["groupings"][idx_tup]["identity"] == "triumvirate":
			print(f"Inside triumvirate", flush=True)
			to_probe = ["positive", "negative"]
			BINODALS["groupings"][idx_tup]["center"]["binodals"] = [np.empty((0,3)), np.empty((0,3))]

			# this portion is to "raise the roof" and define all the good points 
			# start with the central critical point and get all the points where we will find our solutions
			tv_center = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["center"]["idx"] ]["tang_vec"]
			nv_center = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["center"]["idx"] ]["norm_vec"]
			cp_center = P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]]
			sign_uidx = np.sign(np.cross(tv_center, unstable_centers[idx_tup[0]][0:2]-cp_center[0:2]))

			# get the points underneath the tangent of the crit point 
			uphi           = phi[:, 0:2]
			adj_uphi       = (uphi-(cp_center[0:2]+0.01*nv_center))/np.linalg.norm(uphi-(cp_center[0:2]+0.01*nv_center), axis=1).reshape(-1,1)
			uphi           = uphi[np.sign(np.cross(tv_center, adj_uphi)) == sign_uidx]

			# get the hull around the good points
			hull      = ConvexHull(uphi[:, 0:2])
			hull_path = Path(uphi[:, 0:2][hull.vertices])

			for probe in to_probe:
				print(f"In probe = {probe}", flush=True)
				# divide points per the normal now
				crit_probe     = P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]]
				nv_probe       = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"]
				adj_uphi       = (uphi[:,0:2]-crit_probe[0:2])/np.linalg.norm(uphi[:,0:2]-crit_probe[0:2], axis=1).reshape(-1,1)
				clock          = np.sign(np.cross(nv_probe, adj_uphi))
				phi_positive   = uphi[clock>=0]
				phi_negative   = uphi[clock<0 ] 
				
				print(f"Hitting the big calcs...", flush=True)
				results    = P.sym_mu_ps.perform_sweep(phi_positive, phi_negative)
				sol1, sol2 = P.sym_mu_ps.binodal_finder_(results[0], results[1], hull_path)
				sol1, kept = ternary.remove_close_rows(sol1, 1e-15)
				sol2       = sol2[kept]

				neg_sol, pos_sol = clean_and_sort(sol1, sol2, crit_probe, nv_probe)	

				if probe == "positive":
					# in the process of refining, I will only keep points that are on the negative side
					# of positive crit point, but strictly positive for the negative crit point i.e. I will only
					# keep points from from the negative solution that are positive wrt to the central critical point
					nv_center     = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["center"]["idx"] ]["norm_vec"]
					adj_neg_sol   = (neg_sol[:,0:2]-P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2])/np.linalg.norm((neg_sol[:,0:2]-P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(nv_center, adj_neg_sol))
					
					# only keep the points that are positive wrt to the central critical point
					to_keep_1     = clock>=0 
					BINODALS["groupings"][idx_tup][probe]["binodals"]       = [neg_sol[to_keep_1], pos_sol[to_keep_1]] 

					# sort the binodals
					dists = np.linalg.norm(BINODALS["groupings"][idx_tup][probe]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]][0:2], axis=1)
					BINODALS["groupings"][idx_tup][probe]["binodals"][0] = BINODALS["groupings"][idx_tup][probe]["binodals"][0][np.argsort(dists)]
					BINODALS["groupings"][idx_tup][probe]["binodals"][1] = BINODALS["groupings"][idx_tup][probe]["binodals"][1][np.argsort(dists)]

					print(f"Hitting add_and_solve...", flush=True)
					amped_neg_sol, amped_pos_sol = add_and_solve(BINODALS["groupings"][idx_tup][probe]["binodals"][0], BINODALS["groupings"][idx_tup][probe]["binodals"][1], P, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"], 50)
					BINODALS["groupings"][idx_tup][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][idx_tup]["positive"]["idx"]], amped_neg_sol))
					BINODALS["groupings"][idx_tup][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][idx_tup]["positive"]["idx"]], amped_pos_sol))
					print(f"Exiting add_and_solve", flush=True)

					# only keep points on the other side of the negative critical point
					nv_negative   = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["negative"]["idx"] ]["norm_vec"]
					adj_sol2      = (neg_sol[:,0:2] - P.crits[BINODALS["groupings"][idx_tup]["negative"]["idx"]][0:2])/np.linalg.norm((neg_sol[:,0:2] - P.crits[BINODALS["groupings"][idx_tup]["negative"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(nv_negative, adj_sol2))
					to_keep_2     = clock <= 0

					BINODALS["groupings"][idx_tup]["center"]["binodals"][0] = np.vstack((BINODALS["groupings"][idx_tup]["center"]["binodals"][0], neg_sol[to_keep_2]))
					BINODALS["groupings"][idx_tup]["center"]["binodals"][1] = np.vstack((BINODALS["groupings"][idx_tup]["center"]["binodals"][1], pos_sol[to_keep_2]))

				elif probe == "negative":
					# in the process of refining, I will only keep points that are on the negative 
					# of positive crit point, but strictly positive for the negative crit point
					nv_center     = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["center"]["idx"] ]["norm_vec"]
					adj_pos_sol   = (pos_sol[:,0:2]-P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2])/np.linalg.norm((pos_sol[:,0:2]-P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(nv_center, adj_pos_sol))
					
					# only keep the points that are negative wrt to the center critical point
					to_keep_1 = clock <= 0
					BINODALS["groupings"][idx_tup][probe]["binodals"] = [neg_sol[to_keep_1], pos_sol[to_keep_1]]

					# sort the binodals
					dists = np.linalg.norm(BINODALS["groupings"][idx_tup][probe]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]][0:2], axis=1)
					BINODALS["groupings"][idx_tup][probe]["binodals"][0] = BINODALS["groupings"][idx_tup][probe]["binodals"][0][np.argsort(dists)]
					BINODALS["groupings"][idx_tup][probe]["binodals"][1] = BINODALS["groupings"][idx_tup][probe]["binodals"][1][np.argsort(dists)]
					
					
					amped_neg_sol, amped_pos_sol = add_and_solve(BINODALS["groupings"][idx_tup][probe]["binodals"][0], BINODALS["groupings"][idx_tup][probe]["binodals"][1], P, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"], 50)
					BINODALS["groupings"][idx_tup][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_neg_sol))
					BINODALS["groupings"][idx_tup][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_pos_sol))
					
					# only keep points on the other side of the positive critical point
					nv_positive   = BINODALS["crit_info"][ BINODALS["groupings"][idx_tup]["positive"]["idx"] ]["norm_vec"]
					adj_sol1      = (pos_sol[:,0:2] - P.crits[BINODALS["groupings"][idx_tup]["positive"]["idx"]][0:2])/np.linalg.norm((pos_sol[:,0:2] - P.crits[BINODALS["groupings"][idx_tup]["positive"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(nv_positive, adj_sol1))
					to_keep_2     = clock >= 0
					
					BINODALS["groupings"][idx_tup]["center"]["binodals"][0] = np.vstack((BINODALS["groupings"][idx_tup]["center"]["binodals"][0], neg_sol[to_keep_2]))
					BINODALS["groupings"][idx_tup]["center"]["binodals"][1] = np.vstack((BINODALS["groupings"][idx_tup]["center"]["binodals"][1], pos_sol[to_keep_2]))
			
			print(f"Populating the center...", flush=True)
			# populate the centers...
			populate_center(BINODALS, idx_tup, 50)

			# now, make sure after populating the center, make sure there are no gaps
			probe = "negative"
			print(BINODALS["groupings"][idx_tup][probe]["idx"])
			amped_neg_sol, amped_pos_sol = add_and_solve (BINODALS["groupings"][idx_tup][probe]["binodals"][0], BINODALS["groupings"][idx_tup][probe]["binodals"][1], P, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"], 50)
			amped_neg_sol, amped_pos_sol = clean_and_sort(amped_neg_sol, amped_pos_sol, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"])

			BINODALS["groupings"][idx_tup][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_neg_sol))
			BINODALS["groupings"][idx_tup][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_pos_sol))

			combined = np.vstack(BINODALS["groupings"][idx_tup][probe]["binodals"])
			hull      = ConvexHull(combined[:,0:2])
			hull_path = Path(combined[:,0:2][hull.vertices])
			BINODALS["hull_info"]["function"].append(hull_path)
			BINODALS["hull_info"]["binodal"].append([BINODALS["groupings"][idx_tup][probe]["binodals"], "two_phase"])
			
			probe = "positive"
			amped_neg_sol, amped_pos_sol = add_and_solve (BINODALS["groupings"][idx_tup][probe]["binodals"][0], BINODALS["groupings"][idx_tup][probe]["binodals"][1], P, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"], 50)
			amped_neg_sol, amped_pos_sol = clean_and_sort(amped_neg_sol, amped_pos_sol, P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][idx_tup][probe]["idx"] ]["norm_vec"])
			BINODALS["groupings"][idx_tup][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_neg_sol))
			BINODALS["groupings"][idx_tup][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][idx_tup][probe]["idx"]], amped_pos_sol))
			
			combined = np.vstack(BINODALS["groupings"][idx_tup][probe]["binodals"])
			hull      = ConvexHull(combined[:,0:2])
			hull_path = Path(combined[:,0:2][hull.vertices])
			BINODALS["hull_info"]["function"].append(hull_path) 
			BINODALS["hull_info"]["binodal"].append([BINODALS["groupings"][idx_tup][probe]["binodals"], "two_phase"])
			
			# clean up the rest 
			BINODALS["groupings"][idx_tup]["center"]["binodals"][0], keep = ternary.remove_close_rows(BINODALS["groupings"][idx_tup]["center"]["binodals"][0], 1e-15)
			BINODALS["groupings"][idx_tup]["center"]["binodals"][1]       = BINODALS["groupings"][idx_tup]["center"]["binodals"][1][keep]

			dists = np.linalg.norm(BINODALS["groupings"][idx_tup]["center"]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][idx_tup]["center"]["idx"]][0:2], axis=1)
			BINODALS["groupings"][idx_tup]["center"]["binodals"][0] = BINODALS["groupings"][idx_tup]["center"]["binodals"][0][np.argsort(dists)]
			BINODALS["groupings"][idx_tup]["center"]["binodals"][1] = BINODALS["groupings"][idx_tup]["center"]["binodals"][1][np.argsort(dists)]

			combined = np.vstack(BINODALS["groupings"][idx_tup]["center"]["binodals"])
			hull      = ConvexHull(combined[:,0:2])
			hull_path = Path(combined[:,0:2][hull.vertices])
			BINODALS["hull_info"]["function"].append(hull_path) 
			BINODALS["hull_info"]["binodal"].append([BINODALS["groupings"][idx_tup]["center"]["binodals"], "two_phase"])

			B = BINODALS["groupings"][idx_tup]["center"]["binodals"]

			t1 = np.mean([BINODALS["groupings"][idx_tup]["negative"]["binodals"][1][-1], BINODALS["groupings"][idx_tup]["positive"]["binodals"][0][-1]], axis=0)
			t2 = np.mean([BINODALS["groupings"][idx_tup]["negative"]["binodals"][0][-1], BINODALS["groupings"][idx_tup]["center"]  ["binodals"][0][0 ]], axis=0)
			t3 = np.mean([BINODALS["groupings"][idx_tup]["positive"]["binodals"][1][-1], BINODALS["groupings"][idx_tup]["center"]  ["binodals"][1][0 ]], axis=0)
			
			roots = fsolve(triangle_finder, [t1[0], t1[1], t2[0], t2[1], t3[0], t3[1]])
			
			closeness = np.isclose(triangle_finder(roots), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
			
			if closeness.all():
				print("Worked!")
			else: 
				print("This ain't a triangle...", flush=True)

			triangular_hull = np.array([[roots[0], roots[1]], [roots[2], roots[3]], [roots[4], roots[5]]])

			hull = ConvexHull(triangular_hull[:,0:2])
			hull_path = Path(triangular_hull[:,0:2][hull.vertices])
			BINODALS["hull_info"]["triangles"].append(hull_path)
			BINODALS["hull_info"]["function"].append(hull_path)
			BINODALS["hull_info"]["binodal"].append([triangular_hull, "three_phase"])
			
			
	# start the island scans
	if len(stable_islands)>1:
		'''
		for i in range(len(stable_islands)):
			for j in range(i+1, len(stable_islands)):
				# scan the islands
				island_checks = P.sym_mu_ps.perform_sweep(stable_islands[i], stable_islands[j])

				# perform the sweep to get the concrete binodal
				arm_1, arm_2 = P.sym_mu_ps.binodal_finder(island_checks[0], island_checks[1], hull_paths_s[i], hull_paths_s[j])

				# define the dictionary 
				BINODALS["island_scans"][(i, j)] = dict()
				BINODALS["island_scans"][(i, j)][i] = dict() 
				BINODALS["island_scans"][(i, j)][j] = dict() 
				arm_1, kept = ternary.remove_close_rows(arm_1, 1e-15)
				if len(kept) > 3:
					print(f"(i,j) = {(i,j)}")
					print(f"kept = {kept}")
					arm_2 = arm_2[kept]
					BINODALS["island_scans"][(i, j)][i]["binodal"] = arm_1
					BINODALS["island_scans"][(i, j)][j]["binodal"] = arm_2
					idx_tup = (i, j)
					extend_binodals(BINODALS, P, idx_tup)

				else:
					BINODALS["island_scans"][(i, j)][i]["binodal"] = np.empty((0,3))
					BINODALS["island_scans"][(i, j)][j]["binodal"] = np.empty((0,3))
					
		keys = list(BINODALS["island_scans"].keys())

		# sort the binodals
		for key in keys:
			print(f"key = {key}", flush=True)
			print(f'arm1 = {BINODALS["island_scans"][key][key[0]]["binodal"].shape[0]}')
			print(f'arm2 = {BINODALS["island_scans"][key][key[1]]["binodal"].shape[0]}') 
			if BINODALS["island_scans"][key][key[0]]["binodal"].shape[0] != 0:

				# clean out the binodals using remove close rows 
				BINODALS["island_scans"][key][key[0]]["binodal"], kept = ternary.remove_close_rows(BINODALS["island_scans"][key][key[0]]["binodal"], 1e-15)
				BINODALS["island_scans"][key][key[1]]["binodal"]       = BINODALS["island_scans"][key][key[1]]["binodal"][kept]

				# time to sort the binodals
				center = np.mean(np.vstack((BINODALS["island_scans"][key][key[0]]["binodal"], BINODALS["island_scans"][key][key[1]]["binodal"])), axis=0)
				angles = np.arctan2(BINODALS["island_scans"][key][key[0]]["binodal"][:,1] - center[1], BINODALS["island_scans"][key][key[0]]["binodal"][:,0] - center[0])
				
				# sort curve 
				BINODALS["island_scans"][key][key[0]]["binodal"] = BINODALS["island_scans"][key][key[0]]["binodal"][np.argsort(angles)] 
				BINODALS["island_scans"][key][key[1]]["binodal"] = BINODALS["island_scans"][key][key[1]]["binodal"][np.argsort(angles)]

				# find the biggest angular difference 
				angles          = angles[np.argsort(angles)]
				max_delta_theta = np.max(angles[1:] - angles[:-1])

				if max_delta_theta > np.pi/10:
					dists        = np.linalg.norm(BINODALS["island_scans"][key][key[0]]["binodal"][1:] - BINODALS["island_scans"][key][key[0]]["binodal"][:-1], axis=1)
					max_dist_idx = np.argmax(dists)
					BINODALS["island_scans"][key][key[0]]["binodal"] = np.vstack((BINODALS["island_scans"][key][key[0]]["binodal"][max_dist_idx+1:], BINODALS["island_scans"][key][key[0]]["binodal"][:max_dist_idx+1]))
					BINODALS["island_scans"][key][key[1]]["binodal"] = np.vstack((BINODALS["island_scans"][key][key[1]]["binodal"][max_dist_idx+1:], BINODALS["island_scans"][key][key[1]]["binodal"][:max_dist_idx+1]))
				
				# the curve has been sorted
				# now, i would like the curve to begin from the point where the volume fraction is zero, or closest to it
				if BINODALS["island_scans"][key][key[0]]["binodal"][0,0] < 5e-2 or BINODALS["island_scans"][key][key[0]]["binodal"][0,1] < 5e-2 or BINODALS["island_scans"][key][key[0]]["binodal"][0,2] < 5e-2:
					pass
				elif BINODALS["island_scans"][key][key[0]]["binodal"][-1,0] < 5e-2 or BINODALS["island_scans"][key][key[0]]["binodal"][-1,1] < 5e-2 or BINODALS["island_scans"][key][key[0]]["binodal"][-1,2] < 5e-2:
					BINODALS["island_scans"] [key][key[0]]["binodal"] = np.flip(BINODALS["island_scans"][key][key[0]]["binodal"], axis=0)
					BINODALS["island_scans"] [key][key[1]]["binodal"] = np.flip(BINODALS["island_scans"][key][key[1]]["binodal"], axis=0)
				else:
					pass

			else:
				continue
		
		f = open(args.fb, 'wb')
		pickle.dump(BINODALS, f)
		f.close()
		
		print(f"Made the intermediate image...", flush=True)
		for key in BINODALS["island_scans"]:
			B1 = BINODALS["island_scans"][key][key[0]]["binodal"]
			B2 = BINODALS["island_scans"][key][key[1]]["binodal"]
			ax.scatter(B1[:,0], 1-B1[:,0]-B1[:,1], B1[:,1], c='black', s=0.5)
			ax.scatter(B2[:,0], 1-B2[:,0]-B2[:,1], B2[:,1], c='white', s=0.5)

		fig.savefig("debug", dpi=1200, bbox_inches="tight")
		
		
		# process the intersections
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
							BINODALS["island_scans"][keys[i]]["intersections"][keys[j]]["boolean"]    = False
							BINODALS["island_scans"][keys[j]]["intersections"][keys[i]]["boolean"]    = False

				elif keys[i][1] in keys[j]:
					BINODALS["island_scans"][keys[i]]["intersections"][keys[j]] = dict()
					BINODALS["island_scans"][keys[j]]["intersections"][keys[i]] = dict()

					idx_in_j = keys[j].index(keys[i][1])
					L1 = LineString(BINODALS["island_scans"][keys[i]][keys[i][1]]["binodal"][:, 0:2])
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
						# print(f"key i = {keys[i][1]}, key j = {keys[j][idx_in_j]}, intersection = {intersection}.", flush=True)
						# ax.scatter(intersection.x, 1-intersection.y-intersection.x, intersection.y, c='k', s=1, zorder=10)
					
				else:
					# no common island, no intersections 
					print(f"No common island, no intersections. Jump past...", flush=True)
					pass
		
		# let's fill up the binodals
		for key in BINODALS["island_scans"]:
			if BINODALS["island_scans"][key][key[0]]["binodal"].shape[0] == 0 or BINODALS["island_scans"][key][key[1]]["binodal"].shape[0] == 0:
				continue
			else:
				populate_between_islands(BINODALS, key)
		'''
		f = open(args.fb, 'wb')
		pickle.dump(BINODALS, f)
		f.close()

		# now, let's clip the binodals.
		keys = list(BINODALS["island_scans"])
		triangle_set  = set()

		# get the sets of points that make the triangles
		for key in BINODALS["island_scans"]:
			kcount = 0
			triangle = [key]
			for k2 in BINODALS["island_scans"][key]["intersections"]:
				if BINODALS["island_scans"][key]["intersections"][k2]["boolean"] and BINODALS["island_scans"][k2]["intersections"][key]["boolean"]:
					triangle.append(k2)

			if len(triangle) == 3:
				if BINODALS["island_scans"][triangle[0]]["intersections"][triangle[1]]["boolean"] \
					and BINODALS["island_scans"][triangle[0]]["intersections"][triangle[2]]["boolean"] \
					and BINODALS["island_scans"][triangle[1]]["intersections"][triangle[2]]["boolean"]:	
					triangle.sort()
					triangle = tuple(triangle)
					triangle_set.add(triangle)
				else:
					pass

		# go to each triangle, and start finding the three points
		for triangle_keys in triangle_set:
			# get the three points 
			triangle_points = list() 
			for i in range(len(triangle_keys)):
				for j in range(i+1, len(triangle_keys)):
					triangle_points.append(BINODALS["island_scans"][triangle_keys[i]]["intersections"][triangle_keys[j]]["poi"][0:2])
			
			root = fsolve(triangle_finder, [triangle_points[0][0], triangle_points[0][1], triangle_points[1][0], triangle_points[1][1], triangle_points[2][0], triangle_points[2][1]])

			print(f"root = {root}", flush=True)
			closeness = np.isclose(triangle_finder(root), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
			print(f"closeness = {closeness}", flush=True)
			if closeness.all():
				print("Worked!")
				triangular_hull = np.array([[root[0], root[1]], [root[2], root[3]], [root[4], root[5]]])
				hull = ConvexHull(triangular_hull[:,0:2])
				hull_path = Path(triangular_hull[:,0:2][hull.vertices])
				BINODALS["hull_info"]["triangles"].append(hull_path)
				BINODALS["hull_info"]["function"].append(hull_path)
				BINODALS["hull_info"]["binodal"].append([triangular_hull, "three_phase"])

			else: 
				print("This ain't a triangle...", flush=True)

		
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

	# get the binodal on a file
	f = open(args.fb, 'wb')
	pickle.dump(BINODALS, f)
	f.close()

	# now, i need to go into each binodal, and delete off the points 
	# whenever there is a tie-line/points inside the triangle 
	
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
			if BINODALS["hull_info"]["binodal"][idx][0][0].shape[0] == 0 or BINODALS["hull_info"]["binodal"][idx][0][1].shape[0] == 0:
				BINODALS["hull_info"]["function"][idx] = None
			else:
				combined  = np.vstack((BINODALS["hull_info"]["binodal"][idx][0][0], BINODALS["hull_info"]["binodal"][idx][0][1]))
				hull      = ConvexHull(combined[:,0:2])
				hull_path = Path(combined[:,0:2][hull.vertices])
				BINODALS["hull_info"]["function"][idx] = hull_path

		else:
			continue
	
	'''
	# make sure there are no overlapping binodals
	for idx in range(len(BINODALS["hull_info"]["binodal"])):
		if BINODALS["hull_info"]["binodal"][idx][-1] == "three_phase" or (BINODALS["hull_info"]["binodal"][idx][0][0].shape[0] == 0 and BINODALS["hull_info"]["binodal"][idx][0][1].shape[0] == 0):
			continue
		to_keep = list()
		for jdx in range(len(BINODALS["hull_info"]["binodal"][idx][0][0])):
			check = True
			line = np.linspace(BINODALS["hull_info"]["binodal"][idx][0][0][jdx][0:2], BINODALS["hull_info"]["binodal"][idx][0][1][jdx][0:2], 1000)
			for kdx in range(idx+1, len(BINODALS["hull_info"]["function"])):
				if BINODALS["hull_info"]["binodal"][kdx][-1] == "three_phase" or (BINODALS["hull_info"]["binodal"][kdx][0][0].shape[0] == 0 and BINODALS["hull_info"]["binodal"][kdx][0][1].shape[0] == 0):
					continue
				if (BINODALS["hull_info"]["function"][kdx].contains_points(line)).any():
					check = False
					break 
				else:
					continue
			if check:
				to_keep.append(jdx)
		BINODALS["hull_info"]["binodal"][idx][0][0] = BINODALS["hull_info"]["binodal"][idx][0][0][to_keep]
		BINODALS["hull_info"]["binodal"][idx][0][1] = BINODALS["hull_info"]["binodal"][idx][0][1][to_keep]
		combined = np.vstack((BINODALS["hull_info"]["binodal"][idx][0][0], BINODALS["hull_info"]["binodal"][idx][0][1]))
		hull = ConvexHull(combined[:,0:2])
		hull_path = Path(combined[:,0:2][hull.vertices])
		BINODALS["hull_info"]["function"][idx] = hull_path
	'''


	# I have all the binodals on me now. 
	# get the binodal on a file
	f = open(args.fb, 'wb')
	pickle.dump(BINODALS, f)
	f.close()

	# I will start creating the hulls now. 
	# first, we hit the groups. 

	for idx, H in enumerate(BINODALS["hull_info"]["binodal"]):
		print(f"arm1 = {H[0][0][-10:]}")
		print(f"arm2 = {H[0][1][-10:]}")
		if H[0][0].shape[0] == 0 and H[0][1].shape[0] == 0:
			continue
		elif H[-1] == "two_phase":
			ax.scatter(H[0][0][:,0], 1-H[0][0][:,0]-H[0][0][:,1], H[0][0][:,1], s=0.5, c='black')
			ax.scatter(H[0][1][:,0], 1-H[0][1][:,0]-H[0][1][:,1], H[0][1][:,1], s=0.5, c='white')
		elif H[-1] == "three_phase":
			ax.plot(np.hstack([H[0][:,0],H[0][0,0]]),\
			np.hstack([1-H[0][:,0]-H[0][:,1], 1-H[0][0,0]-H[0][0,1]]), np.hstack([H[0][:,1], H[0][0,1]]), c='slategray', lw=1)

	# create the image
	print("Making image...", end=' ')
	
	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		fig.savefig (img_name, dpi=1200, bbox_inches="tight")
	
	elif args.img[-1] == '/':
		fig.savefig(args.img+f"bin_tern-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	
	else:
		fig.savefig (f"bin_tern-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	
	print(f"done!", flush=True)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

