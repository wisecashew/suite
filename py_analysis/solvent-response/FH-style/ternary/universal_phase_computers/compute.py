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
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float, action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float, action='store', help='specific volume of polymer.')
parser.add_argument('--pkl',        dest='pkl',   action='store', type=str,   help='extract information from the pickle file.')
parser.add_argument('--plot-edges', dest='pe',    action='store_true', help='plot the edges of the spinodal.')
parser.add_argument('--plot-crits', dest='pc',    action='store_true', help='plot the critical points.')
parser.add_argument('--no-rtw',     dest='nrtw',  action='store_true',  default=False, help="Dont print out the runtime warning.")
parser.add_argument('--img',        dest='img',   action='store',      type=str, default="None", help="name of image to be created. (default: blastradius).")
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

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
	# crits = np.array([[0.41666667, 0.41666667, 0.16666667],[0.16666667, 0.41666667, 0.41666666], [0.41666667, 0.16666667, 0.41666667]])
	# P.spinodal.crits = crits
	P.spinodal.obtain_crits()
	P.crits = P.spinodal.crits
	print(f"crits = {P.crits}", flush=True)
	print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	print("Plotted out second sweep!", flush=True)
	print(f"Now, we begin island hunts.", flush=True)

	# extract island information
	f = open(args.pkl, 'rb')
	islands = pickle.load(f)
	f.close()

	print(f"Number of islands is {len(islands)}...")

	if len(islands) == 1:
		print(f"Tangent tracing ought to do the job.")
		P.tangent_tracing(ax)

	else:
		# P.tangent_tracing(ax)
		# P.island_hunt()
		# let's start processing the islands. 
		hull_paths = []
		for idx, island in enumerate(islands):
			islands[idx] = np.array(np.vstack([0.001+(0.999-0.001)/999*island[:,1], 0.001+(0.999-0.001)/999*island[:,0]])).T 
			hull = ConvexHull(islands[idx])
			hull_paths.append(Path(islands[idx][hull.vertices]))
		
		# find the island in the middle
		# the island whose center is closest to [1/3, 1/3] is the central island. 
		c_dists = []
		center  = np.array([1/3,1/3])
		for i in range(len(islands)):
			island_c = np.mean(islands[i], axis=0)
			c_dists.append(np.linalg.norm(center-island_c))
		
		c_idx = np.argmin(c_dists)

		binodals  = dict() 
		central = []
		central_curves = []


		for i in range(len(islands)):
			for j in range(i+1, len(islands)):
				results = P.sym_mu_ps.perform_sweep(islands[i], islands[j])
				# perform the sweep 
				sol1, sol2 = P.sym_mu_ps.binodal_finder(results[0], results[1], hull_paths[i], hull_paths[j])

				binodals[(i,j)] = [sol1, sol2]
				if c_idx in (i,j):
					central.append((i,j))
					if i == c_idx:
						central_curves.append(sol1)
					else:
						central_curves.append(sol2)

				# ax.scatter(sol1[:,0], 1-sol1[:,0]-sol1[:,1], sol1[:,1], c="crimson", s=0.5)
				# ax.scatter(sol2[:,0], 1-sol2[:,0]-sol2[:,1], sol2[:,1], c="black", s=0.5)


		# start curating curves
		# sort curves using theta 
		# to do so, find the center of the curves
		center = np.mean(np.vstack((central_curves[0], central_curves[1], central_curves[2])), axis=0)

		# now, start sorting wrt theta
		linestrings  = []
		for c_idx, c_curve in enumerate(central_curves):
			# divide the curves into positive and negative cross product
			direction  = (c_curve - center)/np.linalg.norm(c_curve - center, axis=1)[:, np.newaxis]
			central_axis  = (c_curve[0]-center)/np.linalg.norm(c_curve[0]-center)
			cross_prod    = np.sign(np.cross(direction, central_axis))
			positive_half = (cross_prod == 1)
			negative_half = (cross_prod == -1)

			smack_zero    = c_curve[np.logical_and(cross_prod == 0, np.dot(direction, central_axis)==1) ]
			smack_180     = c_curve[np.logical_and(cross_prod == 0, np.dot(direction, central_axis)==-1)]

			# deal with the positive half
			normalized_pos_half     = (c_curve[positive_half]-center)/np.linalg.norm(c_curve[positive_half]-center)
			pos_theta               = np.arccos(np.dot(normalized_pos_half, central_axis))
			sorted_indices          = np.flip(np.argsort(pos_theta))
			pos_half                = c_curve[positive_half][sorted_indices] 

			# deal with the negative half
			normalized_neg_half     = (c_curve[negative_half]-center)/np.linalg.norm(c_curve[negative_half]-center)
			neg_theta               = np.arccos(np.dot(normalized_neg_half, central_axis))
			sorted_indices          = np.argsort(neg_theta)
			neg_half                = c_curve[negative_half][sorted_indices] 

			# put it all together 
			c_curve = np.vstack((pos_half, smack_zero, pos_neg, smack_180))
			central_curves[c_idx] = c_curve 

			linestrings.append(LineString(c_curve))

		ls_len = len(linestrings)

		for i in range(ls_len):
			for j in range(i+1, ls_len):
				P_intersect = linestrings[i].intersection(linestrings[j])
				# figure out which curve comes first, per theta 
				d1 = (central_curves[i][0] - center)/np.linalg.norm(central_curves[i][0] - center)
				d2 = (central_curves[j][0] - center)/np.linalg.norm(central_curves[j][0] - center)
				t =  np.cross(d1, d2) 
				if t > 0:
					# d2 comes before d1
					intrsxn = linestrings[i].intersection(linestrings[j])
					p       = np.array([intrsxn.x, intrsxn.y])
					# kill everything before p for d2, everything arter p for d1
					axis = (p-center)/np.linalg.norm(p-center)
					direction = (central_curves[i]-center)/np.linalg.norm(central_curves[i]-center, axis=1)
					kept = direction[np.cross(axis, direction)>=0]
					central_curves[i] = central_curves[i][kept]

					direction = (central_curves[j]-center)/np.linalg.norm(central_curves[j]-center, axis=1)
					kept = direction[np.cross(axis, direction)<=0]
					central_curves[j] = central_curves[j][kept]
				
				else:
					# d2 comes before d1
					intrsxn = linestrings[i].intersection(linestrings[j])
					p       = np.array([intrsxn.x, intrsxn.y])
					# kill everything before p for d2, everything arter p for d1
					axis = (p-center)/np.linalg.norm(p-center)
					direction = (central_curves[j]-center)/np.linalg.norm(central_curves[j]-center, axis=1)
					kept = direction[np.cross(axis, direction)>=0]
					central_curves[j] = central_curves[j][kept]

					direction = (central_curves[i]-center)/np.linalg.norm(central_curves[i]-center, axis=1)
					kept = direction[np.cross(axis, direction)<=0]
					central_curves[i] = central_curves[i][kept]

		for curve in central_curves:
			ax.scatter(curve[:,0], 1-curve[:,0]-curve[:,1], curve[:,1], c='gold')
			

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

