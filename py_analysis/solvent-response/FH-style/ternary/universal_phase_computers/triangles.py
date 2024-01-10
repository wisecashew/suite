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
parser.add_argument('--island-stable-pkl',   metavar='SPKL',    dest='spkl',          type=str,   action='store', help='extract information about the stable islands from the pickle file (default: None).',   default=None)
parser.add_argument('--island-unstable-pkl', metavar='UPKL',    dest='upkl',          type=str,   action='store', help='extract information about the unstable islands from the pickle file (default: None).', default=None)
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
			print(f"eq1 = {eq1}, eq2 = {eq2}, eq3 = {eq3}", flush=True)
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

				elif np.lingalg.norm(p1[0:2]-p2[0:2]) < 1e-6:
					continue

				else:
					if np.sign(np.cross(central_axis, p1[0:2]-center)) == np.sign(np.cross(central_axis, p2[0:2]-center)):
						continue
					elif np.cross(central_axis, p1[0:2]-center) >=0:
						b1 = np.vstack((b1, p1))
						b2 = np.bstack((b2, p2))
					else:
						b1 = np.vstack((b1, p2))
						b2 = np.bstack((b2, p1))
					break

	return b1, b2
				
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

	BINODALS = dict()
	BINODALS["hull_info"] = dict()
	BINODALS["groupings"] = dict()
	BINODALS["crit_info"] = dict()

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
		if len(BINODALS["groupings"][uidx]["raw_list"]) == 1:
			BINODALS["groupings"][uidx]["center"] = dict()
			BINODALS["groupings"][uidx]["center"]["idx"] = BINODALS["groupings"][uidx]["raw_list"][0]

		elif len(BINODALS["groupings"][uidx]["raw_list"]) == 3:
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

	# end of setup. Now to move into specifics. 

	if len(stable_islands) == 1 and len(unstable_islands) == 3 and len(P.crits) == 3:
		print(f"Tangent tracing ought to do the job.", flush=True)
		P.tangent_tracing(BINODALS)
		for idx, c in enumerate(P.crits):
			for test_idx in BINODALS["groupings"]:
				if idx in BINODALS["groupings"][test_idx]["raw_list"]:
					uidx = test_idx
					break

			tang_slope   = tangent.tangent2(P.vs, P.vc, P.vp, c[0], c[1], P.chi_pc, P.chi_ps, P.chi_sc, P.spinodal.root_up_s, P.spinodal.root_lo_s)
			normal_slope = -1/tang_slope 
			norm_vec     = np.array([1, normal_slope])/np.sqrt(1+normal_slope**2)
			test_point   = c[0:2] + 0.01 * norm_vec
			if ternary.stab_crit(test_point[0], test_point[1], P.vs, P.vc, P.vp, P.chi_ps, P.chi_pc, P.chi_sc) < 0:
				pass
			else:
				norm_vec = -norm_vec
			neg_sol, pos_sol = clean_and_sort(BINODALS[idx][0], BINODALS[idx][1], c, norm_vec)
			BINODALS["groupings"][uidx]["center"]["binodals"][0] = neg_sol
			BINODALS["groupings"][uidx]["center"]["binodals"][1] = pos_sol
			ax.scatter(neg_sol[:,0], 1-neg_sol[:,0]-neg_sol[:,1], neg_sol[:,1], c='black', s=0.5)
			ax.scatter(pos_sol[:,0], 1-pos_sol[:,0]-pos_sol[:,1], pos_sol[:,1], c='white', s=0.5)
		
		for idx, c in enumerate(P.crits):
			for test_idx in BINODALS["groupings"]:
				if idx in BINODALS["groupings"][test_idx]["raw_list"]:
					uidx = test_idx
					break
			BINODALS["hull_info"]["numerics"] = list()
			BINODALS["hull_info"]["hull"]     = list()
			BINODALS["hull_info"]["splitter"] = list()
			BINODALS["hull_info"]["numerics"].append(BINODALS["groupings"][uidx]["center"]["binodals"])
			stacked = np.vstack((BINODALS["hull_info"]["numerics"][0][:,0:2], BINODALS["hull_info"]["numerics"][1][:,0:2]))
			hull    = ConvexHull(stacked)
			hull    = Path(stacked[hull.vertices])
			BINODALS["hull_info"]["hull"].append(hull)
			def find_split(P):
				return find_closest_line_segment(BINODALS["hull_info"]["numerics"][0][:,0:2], BINODALS["hull_info"]["numerics"][1][:,0:2], P)
			BINODALS["hull_info"]["splitter"].append(find_split)




	elif len(stable_islands) == 2:
		pass 

	elif len(stable_islands) == 4 and len(unstable_islands) == 1 and len(P.crits) == 3:

		# transform the islands and get the hulls
		print(f"Getting the hulls for the islands...", flush=True, end=' ')
		hull_paths = transform_islands(islands) 
		print(f"done!", flush=True)

		# find the island in the middle
		# the island whose center is closest to [1/3, 1/3] is the central island. 

		try:
			print(f"Extracting binodal information from a pickled file...", flush=True, end=' ')
			f = open(args.bpkl, 'rb')
			BINODALS = pickle.load(f)
			f.close()

			print("done!", flush=True)

		except:
			print(f"No file was found, generating a new binodal pkl file.", flush=True)
			print(f"Computing the binodals and other stuff...", flush=True)
			f = open(args.bpkl, 'wb')
			BINODALS              = dict() 

			for i in range(len(islands)):
				for j in range(i+1, len(islands)):
					results = P.sym_mu_ps.perform_sweep(islands[i], islands[j])

					# perform the sweep 
					sol1, sol2 = P.sym_mu_ps.binodal_finder(results[0], results[1], hull_paths[i], hull_paths[j])

					BINODALS[(i, j)]    = dict()
					BINODALS[(i, j)][i] = dict()
					BINODALS[(i, j)][j] = dict()

					BINODALS[(i, j)][i]["points"] = sol1 
					BINODALS[(i, j)][j]["points"] = sol2
					
			# b_info = [binodals, keys_to_center, curves_inside_center, b_linestrings]
			pickle.dump(BINODALS, f)
			f.close()
			print("done!", flush=True)

		# sort out the binodals
		
		print(f"Sorting out all the binodals...", end=' ', flush=True)
		keys = list(BINODALS.keys())
		for key in keys:

			BINODALS[key][key[0]]["points"], kept = ternary.remove_close_rows(BINODALS[key][key[0]]["points"], 1e-6)
			BINODALS[key][key[1]]["points"] = BINODALS[key][key[1]]["points"][kept]

			center = np.mean(np.vstack((BINODALS[key][key[0]]["points"], BINODALS[key][key[1]]["points"])), axis=0)
			angles = np.arctan2(BINODALS[key][key[0]]["points"][:,1] - center[1], BINODALS[key][key[0]]["points"][:,0] - center[0])
			
			# sort curve 
			BINODALS[key][key[0]]["points"] = BINODALS[key][key[0]]["points"][np.argsort(angles)] 
			BINODALS[key][key[1]]["points"] = BINODALS[key][key[1]]["points"][np.argsort(angles)]

			# find the biggest angular difference 
			angles = angles[np.argsort(angles)]
			max_delta_theta = np.max(angles[1:] - angles[:-1])
			# print(f"For key = {key}, max_delta_theta = {max_delta_theta}", flush=True)
			

			if max_delta_theta > np.pi/2:
				dists = np.linalg.norm(BINODALS[key][key[0]]["points"][1:] - BINODALS[key][key[0]]["points"][:-1], axis=1)
				max_dist_idx = np.argmax(dists)
				BINODALS[key][key[0]]["points"] = np.vstack((BINODALS[key][key[0]]["points"][max_dist_idx+1:], BINODALS[key][key[0]]["points"][:max_dist_idx+1]))
				BINODALS[key][key[1]]["points"] = np.vstack((BINODALS[key][key[1]]["points"][max_dist_idx+1:], BINODALS[key][key[1]]["points"][:max_dist_idx+1]))

			# the curve has been sorted. 
			# now, i would like the curve to begin from the point where the volume fraction is zero, or close to it. 
			if BINODALS[key][key[0]]["points"][0,0] < 5e-2 or BINODALS[key][key[0]]["points"][0,1] < 5e-2 or BINODALS[key][key[0]]["points"][0,2] < 5e-2:
				pass
			elif BINODALS[key][key[0]]["points"][-1,0] < 5e-2 or BINODALS[key][key[0]]["points"][-1,1] < 5e-2 or BINODALS[key][key[0]]["points"][-1,2] < 5e-2:
				BINODALS[key][key[0]]["points"] = np.flip(BINODALS[key][key[0]]["points"], axis=0)
				BINODALS[key][key[1]]["points"] = np.flip(BINODALS[key][key[1]]["points"], axis=0)
			else:
				pass
				
		print(f"done!", flush=True)
		
		# now that I have all the binodals, i am going to start finding the intersections for ALL of them. 
		# take two keys, see if they have a common element 
		
		for i in range(len(keys)):
			BINODALS[keys[i]][keys[i][0]]["intersections"] = dict() 
			BINODALS[keys[i]][keys[i][1]]["intersections"] = dict()
			BINODALS[keys[i]]["intersections"] = dict() 
			
			for j in range(len(keys)):
				if i == j:
					continue
				# print(f"i = {keys[i]}, j ={keys[j]}", flush=True)
				if keys[i][0] in keys[j]:
					idx_in_j = keys[j].index(keys[i][0])
					L1 = LineString(BINODALS[keys[i]][keys[i][0]]["points"][:, 0:2])
					L2 = LineString(BINODALS[keys[j]][keys[j][idx_in_j]]["points"][:, 0:2])
					intersection = L1.intersection(L2)
					if intersection.is_empty:
						print(f"No intersection.", flush=True)
					else:
						BINODALS[keys[i]][keys[i][0]]["intersections"][keys[j]] = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						BINODALS[keys[i]]["intersections"][keys[j]] = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						# print(f"key i = {keys[i][0]}, key j = {keys[j][idx_in_j]}, intersection = {intersection}.", flush=True)

				elif keys[i][1] in keys[j]:
					idx_in_j = keys[j].index(keys[i][1])
					L1 = LineString(BINODALS[keys[i]][keys[i][1]]["points"][:, 0:2])
					L2 = LineString(BINODALS[keys[j]][keys[j][idx_in_j]]["points"][:, 0:2])
					intersection = L1.intersection(L2)
					if intersection.is_empty:
						print(f"No intersection.", flush=True)
					else:
						# print(f"intersection = {intersection}")
						BINODALS[keys[i]][keys[i][1]]["intersections"][keys[j]] = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						BINODALS[keys[i]]["intersections"][keys[j]]             = np.array([intersection.x, intersection.y, 1-intersection.x-intersection.y])
						# print(f"key i = {keys[i][1]}, key j = {keys[j][idx_in_j]}, intersection = {intersection}.", flush=True)
						# ax.scatter(intersection.x, 1-intersection.y-intersection.x, intersection.y, c='k', s=1, zorder=10)
					
				else:
					# no common island, no intersections 
					print(f"Jump past...", flush=True)
					pass
		
		# get the triangles/triads
		BINODALS["triad_curves"] = list()
		BINODALS["triad_points"] = list()
		BINODALS["triad_hulls" ] = list()

		for key in BINODALS:
			if isinstance(key, str):
				continue
			break_out = False 
			# print(f"key = {key}")
			mu_equality = list()
			mu_equality.append(key)

			# now, check if all the binodals that intersect with the binodal of interest intersect with one another
			for i1, e1 in enumerate(BINODALS[key]["intersections"]):
				if break_out:
					break
				for i2, e2 in enumerate(BINODALS[key]["intersections"]):
					# print(f"\te1 = {e1}, e2 = {e2}")
					if i1 >= i2:
						# print("\tMoving on...")
						continue
					elif e2 in BINODALS[e1]["intersections"]:
						mu_equality.append(e1)
						mu_equality.append(e2)
						continue
					else:
						break_out = True
						break
			if not break_out:
				BINODALS["triad_curves"].append(tuple(mu_equality))

		print(f"curves that lead to triads = {BINODALS['triad_curves']}", flush=True)
		for triad in BINODALS["triad_curves"]:
			points = list()
			for i1, k1 in enumerate(triad):
				for i2, k2 in enumerate(triad):
					if i1 >= i2:
						continue
					if k2 in BINODALS[k1]["intersections"]:
						points.append(BINODALS[k1]["intersections"][k2])
			BINODALS["triad_points"].append(points)
				
		print(f"triads = {BINODALS['triad_points']}", flush=True)
		# now that I have all the points, it is time to start cutting up the curves
		for key in BINODALS:

			if isinstance(key, str):
				continue

			# this is to cut up the curves on the sides 
			if len(BINODALS[key]["intersections"]) == 2:
				
				key_list1 = list(BINODALS[key][key[0]]["intersections"].keys())
				key_list2 = list(BINODALS[key][key[1]]["intersections"].keys())

				# print(f"The binodal = {BINODALS[key][key[0]]['points']}")
				# print(f"The point of intersection = {BINODALS[key][key[0]]['intersections'][key_list1[0]][:2]}")

				trunc_index1 = truncate_curve_at_point(BINODALS[key][key[0]]["points"][:,0:2], BINODALS[key][key[0]]["intersections"][key_list1[0]][:2])
				# print(f"trunc_index1 = {trunc_index1}", flush=True)
				# print(f"The cut binodal = {BINODALS[key][key[0]]['points'][:trunc_index1]}")
				BINODALS[key][key[0]]["points"] = np.vstack((BINODALS[key][key[0]]["points"][:trunc_index1], BINODALS[key][key[0]]["intersections"][key_list1[0]]))

				trunc_index2 = truncate_curve_at_point(BINODALS[key][key[1]]["points"][:,0:2], BINODALS[key][key[1]]["intersections"][key_list2[0]][:2])
				# print(f"trunc_index1 = {trunc_index2}", flush=True)
				# print(f"The cut binodal = {BINODALS[key][key[1]]['points'][:trunc_index2]}")
				BINODALS[key][key[1]]["points"] = np.vstack((BINODALS[key][key[1]]["points"][:trunc_index2], BINODALS[key][key[1]]["intersections"][key_list2[0]]))
				if trunc_index1 == trunc_index2:
					print("truncation indices match.")
				else:
					print("truncation indices dont match. Something is wrong.")
			
			elif len(BINODALS[key]["intersections"]) == 4:
				key_list = list(BINODALS[key][key[0]]["intersections"].keys())
				t_ind = list() 
				points = list() 
				for k in key_list:
					points.append(BINODALS[key][key[0]]["intersections"][k])
					t_ind.append(truncate_curve_at_point(BINODALS[key][key[0]]["points"][:,0:2], BINODALS[key][key[0]]["intersections"][k][:2]))
				points = np.array(points)
				t_ind = np.array(t_ind)
				points = points[np.argsort(t_ind)]
				t_ind  = np.sort(t_ind)
				# truncate the central curves
				BINODALS[key][key[0]]["points"] = np.vstack((points[0], BINODALS[key][key[0]]["points"][t_ind[0]:t_ind[1]], points[1]))

				key_list = list(BINODALS[key][key[1]]["intersections"].keys())
				t_ind = list() 
				points = list() 
				for k in key_list:
					points.append(BINODALS[key][key[1]]["intersections"][k])
					t_ind.append(truncate_curve_at_point(BINODALS[key][key[1]]["points"][:,0:2], BINODALS[key][key[1]]["intersections"][k][:2]))
				points = np.array(points)
				t_ind = np.array(t_ind)
				points = points[np.argsort(t_ind)]
				t_ind  = np.sort(t_ind)
				# truncate the central curves
				BINODALS[key][key[1]]["points"] = np.vstack((points[0], BINODALS[key][key[1]]["points"][t_ind[0]:t_ind[1]], points[1]))
		
		# now make every binodal curve and triad into a convex hull
		for key in BINODALS:
			if isinstance(key, str):
				continue
			# print(f"key = {key}")
			# print(np.vstack((BINODALS[key][key[0]]["points"][:,0:2], BINODALS[key][key[1]]["points"][:,0:2])))
			BINODALS[key]["hull"] = ConvexHull(np.vstack((BINODALS[key][key[0]]["points"][:,0:2], BINODALS[key][key[1]]["points"][:,0:2])))

		# got the hull of the triad that allows for three part separation
		for triad in BINODALS["triad_points"]:
			# print(f"triad = {triad}", flush=True)
			BINODALS["triad_hulls"].append(ConvexHull(np.array(triad)[:,0:2]))
		
		
		print(f"keys = {BINODALS.keys()}")
		
		cols = ["bisque", "darkorange", "burlywood"]
		for tidx, triad in enumerate(BINODALS["triad_points"]):
			ax.scatter(triad[0][0], triad[0][2], triad[0][1], marker='o', c=cols[tidx], s=2, zorder=10)
			ax.scatter(triad[1][0], triad[1][2], triad[1][1], marker='o', c=cols[tidx], s=2, zorder=10)
			ax.scatter(triad[2][0], triad[2][2], triad[2][1], marker='o', c=cols[tidx], s=2, zorder=10)

		cols = ["white", "pink", "gold", "grey", "skyblue", "lavender"]
		for bidx, bkey in enumerate(BINODALS):
			if isinstance(bkey, str):
				continue
			for idx in bkey:
				# print(f"BINODALS[{bkey}][{idx}][\"points\"][0:10] = {BINODALS[bkey][idx]['points'][0:10]}", flush=True)
				ax.plot(BINODALS[bkey][idx]["points"][:,0], 1-BINODALS[bkey][idx]["points"][:,0]-BINODALS[bkey][idx]["points"][:,1], BINODALS[bkey][idx]["points"][:,1], lw=0.5, c=cols[bidx%len(cols)])
				# ax.scatter(BINODALS[bkey][idx]["points"][:,0], 1-BINODALS[bkey][idx]["points"][:,0]-BINODALS[bkey][idx]["points"][:,1], BINODALS[bkey][idx]["points"][:,1], s=0.5, c=cols[bidx%len(cols)])
		 		
	f = open(args.fb, 'wb')
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

