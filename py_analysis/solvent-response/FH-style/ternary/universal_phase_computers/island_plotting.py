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
parser.add_argument('--search-density',      metavar='SD',      dest='sd',            type=int,   action='store', help='density of points sampled for stability plot (default: 500).',                          default=500 )
parser.add_argument('--island-stable-pkl',   metavar='SPKL',    dest='spkl',          type=str,   action='store', help='extract information about the stable islands from the pickle file (default: None).',    default=None)
parser.add_argument('--island-unstable-pkl', metavar='UPKL',    dest='upkl',          type=str,   action='store', help='extract information about the unstable islands from the pickle file (default: None).',  default=None)
parser.add_argument('--crit-pkl',            metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                                default=None)
parser.add_argument('--binodal-pkl',         metavar='bpkl',    dest='bpkl',          type=str,   action='store', help='enter name of file with all the information about the binodals (default: None).',       default=None)
parser.add_argument('--plot-edges',          dest='pe',         action='store_true',  help='plot the edges of the spinodal.')
parser.add_argument('--plot-crits',          dest='pc',         action='store_true',  help='plot the critical points.'      )
parser.add_argument('--no-rtw',              dest='nrtw',       action='store_true',  default=False, help="Dont print out the runtime warning."         )
parser.add_argument('--img',                 dest='img',        action='store',       type=str,      default="None", help="name of image to be created.")
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
	# print(f"crits = {P.crits}", flush=True)
	# print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	# P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)
	
	f = open(args.spkl, 'rb')
	stable_islands = pickle.load(f)
	f.close()

	cols_s = ["paleturquoise", "cyan", "darkturquoise", "steelblue"]
	stable_hull_paths = transform_islands(stable_islands) 

	for idx, si in enumerate(stable_islands):
		ax.scatter(si[:,0], 1-si[:,0]-si[:,1], si[:,1], c=cols_s[idx%len(cols_s)], s=1, zorder=5)
	
	print(f"Number of stable islands is {len(stable_islands)}...")

	f = open(args.upkl, 'rb')
	unstable_islands = pickle.load(f)
	print(f"Number of unstable islands is {len(unstable_islands)}")
	f.close()

	cols_u = ["lightcoral", "indianred", "firebrick", "darkred"]
	unstable_hull_paths = transform_islands(unstable_islands)

	for idx, ui in enumerate(unstable_islands):
		ax.scatter(ui[:,0], 1-ui[:,0]-ui[:,1], ui[:,1], c=cols_u[idx%len(cols_u)], s=1, zorder=5)
	
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
		fig.savefig (f"islands_tern-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	print(f"done!", flush=True)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

