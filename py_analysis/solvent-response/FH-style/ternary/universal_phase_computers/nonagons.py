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
	pos_indices = np.where((signs1== 1) & (signs2==-1))[0]
	neg_indices = np.where((signs1==-1) & (signs2== 1))[0]

	# use indices to separate rows into new arrays
	pos_sol = np.vstack((sol1[pos_indices], sol2[neg_indices]))
	neg_sol = np.vstack((sol2[pos_indices], sol1[neg_indices]))

	dists = np.linalg.norm(pos_sol[:,0:2]-crit[0:2], axis=1)
	pos_sol = pos_sol[np.argsort(dists)]
	neg_sol = neg_sol[np.argsort(dists)]

	return pos_sol, neg_sol

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
	pop_array       = np.insert(array, idx+1, insert_rows[1:-1], axis=0)
	return pop_array

#########################################

def solve_within(a1, a2, P, center, central_axis):

	b1 = np.empty((0,3))
	b2 = np.empty((0,3))

	for idx, pt in enumerate(a1):
		def mu_equations(phi):
			eq1 = P.sym_mu_ps.delta_mu_s(pt[0], phi[0], phi[1], phi[2]) # /(np.linalg.norm(np.array([pt[0], phi[0]])-np.array([phi[1], phi[2]])))
			eq2 = P.sym_mu_ps.delta_mu_p(pt[0], phi[0], phi[1], phi[2]) # /(np.linalg.norm(np.array([pt[0], phi[0]])-np.array([phi[1], phi[2]])))
			eq3 = P.sym_mu_ps.delta_mu_c(pt[0], phi[0], phi[1], phi[2]) # /(np.linalg.norm(np.array([pt[0], phi[0]])-np.array([phi[1], phi[2]])))
			# print(f"eq1 = {eq1}, eq2 = {eq2}, eq3 = {eq3}", flush=True)
			return [eq1, eq2, eq3]

		for tidx, tpt in enumerate(a2):
			root = fsolve(mu_equations, [pt[1], tpt[0], tpt[1]])
			if (np.abs(np.array(mu_equations(root)))>1e-6).any():
				continue
			else:
				p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
				p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

				if np.isnan(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isnan(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc)):
					continue

				elif np.linalg.norm(p1[0:2]-p2[0:2]) < 1e-6:
					continue

				else:
					if np.sign(np.cross(central_axis, p1[0:2]-center[0:2])) == np.sign(np.cross(central_axis, p2[0:2]-center[0:2])):
						continue
					elif np.cross(central_axis, p1[0:2]-center[0:2]) >=0:
						b1 = np.vstack((b1, p1))
						b2 = np.vstack((b2, p2))
					else:
						b1 = np.vstack((b1, p2))
						b2 = np.vstack((b2, p1))
					break

	return b1, b2
				
#########################################

def add_and_solve(arm1, arm2, P, center, central_axis, M):
	# find point with greatest gap
	dist1     = np.linalg.norm(np.diff(arm1, axis=0), axis=1)
	max_dist1 = np.argmax(dist1)
	dist2     = np.linalg.norm(np.diff(arm2, axis=0), axis=1)
	max_dist2 = np.argmax(dist2)

	print(f"dist1 = {dist1[max_dist1]}, dist2 = {dist2[max_dist2]}", flush=True)
	iter = 0
	while dist1[max_dist1] > 0.0005 or dist2[max_dist2] > 0.0005:
		
		if dist1[max_dist1] > dist2[max_dist2] and dist1[max_dist1] > 0.0005:
			new_arm1 = add_rows(arm1, M, max_dist1)
			new_arm2 = add_rows(arm2, M, max_dist1)

			b1, b2   = solve_within(new_arm1[max_dist1+1:max_dist1+1+M], new_arm2[max_dist1+1:max_dist1+1+M], P, center, central_axis)
			arm1 = np.vstack((arm1, b1))
			arm2 = np.vstack((arm2, b2))
			arm1, arm2 = clean_and_sort(arm1, arm2, center, central_axis)

		elif dist2[max_dist1] > dist1[max_dist2] and dist2[max_dist2] > 0.0005:
			new_arm1 = add_rows(arm1, M, max_dist2)
			new_arm2 = add_rows(arm2, M, max_dist2)

			b1, b2   = solve_within(new_arm1[max_dist2+1:max_dist2+1+M], new_arm2[max_dist2+1:max_dist2+1+M], P, center, central_axis)
			arm1 = np.vstack((arm1, b1))
			arm2 = np.vstack((arm2, b2))

			arm1, arm2 = clean_and_sort(arm1, arm2, center, central_axis)

		# find point with greatest gap
		dist1     = np.linalg.norm(np.diff(arm1, axis=0), axis=1)
		max_dist1 = np.argmax(dist1)
		dist2     = np.linalg.norm(np.diff(arm2, axis=0), axis=1)
		max_dist2 = np.argmax(dist2)
		arm1, keep = ternary.remove_close_rows(arm1, 1e-12)
		arm2       = arm2[keep]
		print(f"dist1 = {dist1[max_dist1]}, dist2 = {dist2[max_dist2]}", flush=True)
		iter += 1
		if iter > 100:
			break

	return arm1, arm2 




#########################################

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
	
	'''
	f = open('nonagons.pkl', 'rb')
	BINODALS = pickle.load(f)
	f.close() 
	uidx=0

	# print(f'BINODALS[groupings][{uidx}][center][binodals][0] = {BINODALS["groupings"][uidx]["center"]["binodals"][0]}')
	# print(f'BINODALS[groupings][{uidx}][center][binodals][1] = {BINODALS["groupings"][uidx]["center"]["binodals"][1]}')
	# print(f"Critical points along binodal = \n{BINODALS['groupings'][0]['raw_crits'], BINODALS['groupings'][0]['raw_list']}")

	c_cp = P.crits[BINODALS['groupings'][uidx]['center']  ['idx']]
	n_cp = P.crits[BINODALS['groupings'][uidx]['negative']['idx']]
	p_cp = P.crits[BINODALS['groupings'][uidx]['positive']['idx']]
	print(f"Central critical point = {BINODALS['groupings'][uidx]['center']['idx'], c_cp}")
	print(f"Negative critical point = {BINODALS['groupings'][uidx]['negative']['idx'], n_cp}")
	print(f"Positive critical point = {BINODALS['groupings'][uidx]['positive']['idx'], p_cp}")
	# find the triangles
	# choose the center of the positive-negative binodals
	c1 = BINODALS["groupings"][uidx]["positive"]["binodals"][0]
	c2 = BINODALS["groupings"][uidx]["positive"]["binodals"][1]

	c1, c2 = add_and_solve(c1, c2, P, p_cp, BINODALS["crit_info"][ BINODALS['groupings'][uidx]['positive']['idx'] ]["norm_vec"], 50)

	ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
	ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)

	ax.scatter(c_cp[0], c_cp[2], c_cp[1], c='red',s=0.5)
	ax.scatter(p_cp[0], p_cp[2], p_cp[1], c='green', s=0.5)
	ax.scatter(n_cp[0], n_cp[2], n_cp[1], c='yellow', s=0.5)

	fig.savefig(args.img, dpi=1200, bbox_inches='tight')
	exit()
	'''
	BINODALS = dict()
	BINODALS["groupings"] = dict() 
	BINODALS["crit_info"] = dict()
	
	if len(P.crits) == 9:		

		for cidx, crit in enumerate(P.crits):
			BINODALS["crit_info"][cidx] = dict()
			BINODALS["crit_info"][cidx]["tang_slope"] = tangent.tangent2(P.vs, P.vc, P.vp, crit[0], crit[1], P.chi_pc, P.chi_ps, P.chi_sc, P.spinodal.root_up_s, P.spinodal.root_lo_s)
			BINODALS["crit_info"][cidx]["tang_vec"]   = np.array([1, BINODALS["crit_info"][cidx]["tang_slope"]])/np.sqrt(1+BINODALS["crit_info"][cidx]["tang_slope"]**2)
			BINODALS["crit_info"][cidx]["norm_slope"] = -1/BINODALS["crit_info"][cidx]["tang_slope"]
			norm_vec = np.array([1, BINODALS["crit_info"][cidx]["norm_slope"]])/np.sqrt(1+BINODALS["crit_info"][cidx]["norm_slope"]**2)
			test_point = crit[0:2] + 0.05*norm_vec 
			if ternary.stab_crit(test_point[0], test_point[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) > 0:
				BINODALS["crit_info"][cidx]["norm_vec"]   = norm_vec
			else:
				BINODALS["crit_info"][cidx]["norm_vec"]   = -norm_vec

			uidx = np.argmin(np.linalg.norm(unstable_centers - crit[0:2], axis=1))
			if uidx in list(BINODALS["groupings"].keys()):
				BINODALS["groupings"][uidx]["raw_list" ].append(cidx)
				BINODALS["groupings"][uidx]["raw_crits"].append(crit)
			else:
				BINODALS["groupings"][uidx] = dict()
				BINODALS["groupings"][uidx]["binodals"]  = []
				BINODALS["groupings"][uidx]["raw_list"]  = [cidx] 
				BINODALS["groupings"][uidx]["raw_crits"] = [crit]
			
		# time to order the groupings
		for uidx in list(BINODALS["groupings"].keys()):
			BINODALS["groupings"][uidx]["center"  ] = dict()
			BINODALS["groupings"][uidx]["positive"] = dict()
			BINODALS["groupings"][uidx]["negative"] = dict()

			# find the middle point
			com = np.mean(BINODALS["groupings"][uidx]["raw_crits"], axis=0)
			BINODALS["groupings"][uidx]["center"]["idx"] = BINODALS["groupings"][uidx]["raw_list"][np.argmin(np.linalg.norm(np.array(BINODALS["groupings"][uidx]["raw_crits"])[:,0:2]- com[0:2], axis=1))]
			central_axis = BINODALS["crit_info"][BINODALS["groupings"][uidx]["center"]["idx"]]["norm_vec"]

			for c in BINODALS["groupings"][uidx]["raw_list"]:
				if c == BINODALS["groupings"][uidx]["center"]["idx"]:
					continue 
				else:
					deviation = (P.crits[c][0:2] - P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2])/np.linalg.norm(P.crits[c][0:2] - P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2])
					sign = np.sign(np.cross(BINODALS["crit_info"][ BINODALS["groupings"][uidx]["center"]["idx"] ]["norm_vec"][0:2], deviation[0:2]))
					if sign > 0:
						BINODALS["groupings"][uidx]["positive"]["idx"] = c
					else:
						BINODALS["groupings"][uidx]["negative"]["idx"] = c


		# now that I have the critical points grouped and identified as positive, negative, center, I can start constructing binodals
		print(f"keys = {list(BINODALS['groupings'].keys())}")
		for uidx in [0]: # list(BINODALS["groupings"].keys()):
			print(f"uidx = {uidx}", flush=True)
			to_probe = ["positive", "negative"]
			BINODALS["groupings"][uidx]["center"]["binodals"]    = [np.empty((0,3)), np.empty((0,3))] 

			# this portion is to "raise the roof", and define all the good points
			# start with central crit point and get all the points where we will find our solutions
			tangent_vector = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["center"]["idx"] ]["tang_vec"]
			normal_vector  = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["center"]["idx"] ]["norm_vec"]
			central_crit   = P.crits[BINODALS["groupings"][uidx]["center"]["idx"]]
			sign_uidx      = np.sign(np.cross(tangent_vector, unstable_centers[uidx][0:2]-central_crit[0:2]))

			# get the points underneath the tangent of the crit point 
			uphi           = phi[:, 0:2]
			adj_uphi       = (uphi-(central_crit[0:2]+0.05*normal_vector))/np.linalg.norm(uphi-(central_crit[0:2]+0.05*normal_vector), axis=1).reshape(-1,1)
			uphi           = uphi[np.sign(np.cross(tangent_vector, adj_uphi)) == sign_uidx]

			# get the hull around the good points
			hull      = ConvexHull(uphi[:, 0:2])
			hull_path = Path(uphi[:, 0:2][hull.vertices])

			for probe in to_probe:
				
				print(f"Probe in {probe}...", flush=True)

				# divide points per the normal now
				central_crit   = P.crits[BINODALS["groupings"][uidx][probe]["idx"]]
				normal_vector  = BINODALS["crit_info"][ BINODALS["groupings"][uidx][probe]["idx"] ]["norm_vec"]
				adj_uphi       = (uphi[:,0:2]-central_crit[0:2])/np.linalg.norm(uphi[:,0:2]-central_crit[0:2], axis=1).reshape(-1,1)
				clock          = np.sign(np.cross(normal_vector, adj_uphi))
				phi_positive   = uphi[clock>=0]
				phi_negative   = uphi[clock<0] 

				# ax.scatter(phi_positive[:,0], 1-phi_positive[:,0]-phi_positive[:,1], phi_positive[:,1], s=1, c='honeydew')
				# ax.scatter(phi_negative[:,0], 1-phi_negative[:,0]-phi_negative[:,1], phi_negative[:,1], s=1, c='lavender')
				
				print(f"Hitting the big calcs...", flush=True)
				results    = P.sym_mu_ps.perform_sweep(phi_positive, phi_negative)
				sol1, sol2 = P.sym_mu_ps.binodal_finder_(results[0], results[1], hull_path)
				sol1, kept = ternary.remove_close_rows(sol1, 1e-6)
				sol2       = sol2[kept]

				pos_sol, neg_sol = clean_and_sort(sol1, sol2, central_crit, normal_vector)

				# sort by how close it is to the critical point
				dists   = np.linalg.norm(pos_sol[:,0:2]-central_crit[0:2], axis=1)
				pos_sol = pos_sol[np.argsort(dists)]
				neg_sol = neg_sol[np.argsort(dists)]
				
				if probe == "positive":
					# in the process of refining, I will only keep points that are on the negative 
					# of positive crit point, but strictly positive for the negative crit point
					normal_vector = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["center"]["idx"] ]["norm_vec"]
					adj_neg_sol   = (neg_sol[:,0:2]-P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2])/np.linalg.norm((neg_sol[:,0:2]-P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(normal_vector, adj_neg_sol))
					
					# only keep the points that are positive wrt to the center critical point
					to_keep_1       = clock>=0 
					BINODALS["groupings"][uidx][probe]["binodals"]       = [neg_sol[to_keep_1], pos_sol[to_keep_1]] 

					# sort the binodals
					dists = np.linalg.norm(BINODALS["groupings"][uidx][probe]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][uidx][probe]["idx"]][0:2], axis=1)
					BINODALS["groupings"][uidx][probe]["binodals"][0] = BINODALS["groupings"][uidx][probe]["binodals"][0][np.argsort(dists)]
					BINODALS["groupings"][uidx][probe]["binodals"][1] = BINODALS["groupings"][uidx][probe]["binodals"][1][np.argsort(dists)]

					pos_sol, neg_sol = add_and_solve(BINODALS["groupings"][uidx][probe]["binodals"][0], BINODALS["groupings"][uidx][probe]["binodals"][1], P, central_crit, BINODALS["crit_info"][ BINODALS["groupings"][uidx][probe]["idx"] ]["norm_vec"], 50)
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((central_crit, neg_sol))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((central_crit, pos_sol))

					# only keep points on the other side of the negative critical point
					normal_vector = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["negative"]["idx"] ]["norm_vec"]
					adj_sol2      = (neg_sol[:,0:2] - P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]][0:2])/np.linalg.norm((neg_sol[:,0:2] - P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(normal_vector, adj_sol2))
					to_keep_2     = clock <= 0

					BINODALS["groupings"][uidx]["center"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["center"]["binodals"][0], neg_sol[to_keep_2]))
					BINODALS["groupings"][uidx]["center"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["center"]["binodals"][1], pos_sol[to_keep_2]))

					# c1 = BINODALS["groupings"][uidx][probe]["binodals"][0]
					# c2 = BINODALS["groupings"][uidx][probe]["binodals"][1]
					# ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
					# ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)
				
				elif probe == "negative":
					# in the process of refining, I will only keep points that are on the negative 
					# of positive crit point, but strictly positive for the negative crit point
					normal_vector = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["center"]["idx"] ]["norm_vec"]
					adj_pos_sol   = (pos_sol[:,0:2]-P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2])/np.linalg.norm((pos_sol[:,0:2]-P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(normal_vector, adj_pos_sol))
					
					# only keep the points that are negative wrt to the center critical point
					to_keep_1 = clock <= 0
					BINODALS["groupings"][uidx][probe]["binodals"] = [neg_sol[to_keep_1], pos_sol[to_keep_1]]

					# sort the binodals
					dists = np.linalg.norm(BINODALS["groupings"][uidx][probe]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][uidx][probe]["idx"]][0:2], axis=1)
					BINODALS["groupings"][uidx][probe]["binodals"][0] = BINODALS["groupings"][uidx][probe]["binodals"][0][np.argsort(dists)]
					BINODALS["groupings"][uidx][probe]["binodals"][1] = BINODALS["groupings"][uidx][probe]["binodals"][1][np.argsort(dists)]

					# only keep points on the other side of the positive critical point
					normal_vector = BINODALS["crit_info"][ BINODALS["groupings"][uidx]["positive"]["idx"] ]["norm_vec"]
					adj_sol1      = (pos_sol[:,0:2] - P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]][0:2])/np.linalg.norm((pos_sol[:,0:2] - P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]][0:2]), axis=1).reshape(-1,1)
					clock         = np.sign(np.cross(normal_vector, adj_sol1))
					to_keep_2     = clock >= 0
					
					BINODALS["groupings"][uidx]["center"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["center"]["binodals"][0], neg_sol[to_keep_2]))
					BINODALS["groupings"][uidx]["center"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["center"]["binodals"][1], pos_sol[to_keep_2]))

					# c1 = BINODALS["groupings"][uidx][probe]["binodals"][0]
					# c2 = BINODALS["groupings"][uidx][probe]["binodals"][1]
					#ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
					# ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)
			
			# compile the central binodals
			BINODALS["groupings"][uidx]["center"]["binodals"][0], keep = ternary.remove_close_rows(BINODALS["groupings"][uidx]["center"]["binodals"][0], 1e-6)
			BINODALS["groupings"][uidx]["center"]["binodals"][1]       = BINODALS["groupings"][uidx]["center"]["binodals"][1][keep]

			dists = np.linalg.norm(BINODALS["groupings"][uidx]["center"]["binodals"][0][:,0:2]-P.crits[BINODALS["groupings"][uidx]["center"]["idx"]][0:2], axis=1)
			BINODALS["groupings"][uidx]["center"]["binodals"][0] = BINODALS["groupings"][uidx]["center"]["binodals"][0][np.argsort(dists)]
			BINODALS["groupings"][uidx]["center"]["binodals"][1] = BINODALS["groupings"][uidx]["center"]["binodals"][1][np.argsort(dists)]

			c1 = BINODALS["groupings"][uidx]["center"]["binodals"][0]
			c2 = BINODALS["groupings"][uidx]["center"]["binodals"][1]
			ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
			ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)
			
			pos_sol, neg_sol = add_and_solve(BINODALS["groupings"][uidx]["positive"]["binodals"][0], BINODALS["groupings"][uidx]["positive"]["binodals"][1], P, P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][uidx]["positive"]["idx"] ]["norm_vec"], 50)
			BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]], neg_sol))
			BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][uidx]["positive"]["idx"]], pos_sol))

			c1 = BINODALS["groupings"][uidx]["positive"]["binodals"][0]
			c2 = BINODALS["groupings"][uidx]["positive"]["binodals"][1]
			ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
			ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)
			
			pos_sol, neg_sol = add_and_solve(BINODALS["groupings"][uidx]["negative"]["binodals"][0], BINODALS["groupings"][uidx]["negative"]["binodals"][1], P, P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]], BINODALS["crit_info"][ BINODALS["groupings"][uidx]["negative"]["idx"] ]["norm_vec"], 50)
			BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]], neg_sol))
			BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((P.crits[BINODALS["groupings"][uidx]["negative"]["idx"]], pos_sol))

			c1 = BINODALS["groupings"][uidx]["negative"]["binodals"][0]
			c2 = BINODALS["groupings"][uidx]["negative"]["binodals"][1]
			ax.scatter(c1[:,0], 1-c1[:,0]-c1[:,1], c1[:,1], c='black', s=0.5)
			ax.scatter(c2[:,0], 1-c2[:,0]-c2[:,1], c2[:,1], c='white', s=0.5)

			# find the triangles
			# choose the center of the positive-negative binodals
			t1 = (BINODALS["groupings"][uidx]["negative"]["binodals"][1][-1]+BINODALS["groupings"][uidx]["positive"]["binodals"][0][-1])/2

			# choose the center of negative-center binodals
			t2 = (BINODALS["groupings"][uidx]["negative"]["binodals"][0][-1]+BINODALS["groupings"][uidx]["center"]["binodals"][0][0])/2

			# choose the center of positive-center binodals
			t3 = (BINODALS["groupings"][uidx]["positive"]["binodals"][1][-1]+BINODALS["groupings"][uidx]["center"]["binodals"][1][0])/2
			
			t  = np.array([t1, t2, t3, t1])
			ax.plot(t[:,0], 1-t[:,0]-t[:,1], t[:,1], c='green', lw=1)
			

	f = open("nonagons.pkl", 'wb')
	pickle.dump(BINODALS, f)
	f.close()
	
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

