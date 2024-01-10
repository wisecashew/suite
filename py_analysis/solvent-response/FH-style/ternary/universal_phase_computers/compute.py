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
parser.add_argument('--chisc',          metavar='chi_sc',  dest='chi_sc',        type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',          metavar='chi_ps',  dest='chi_ps',        type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',          metavar='chi_pc',  dest='chi_pc',        type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',              metavar='vs',      dest='vs',            type=float, action='store', help='specific volume of solvent.'  )
parser.add_argument('-vc',              metavar='vc',      dest='vc',            type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',              metavar='vp',      dest='vp',            type=float, action='store', help='specific volume of polymer.'  )
parser.add_argument('--search-density', metavar='SD',      dest='sd',            type=int,   action='store', help='density of points sampled for stability plot (default: 500).',                    default=500 )
parser.add_argument('--pkl',            metavar='PKL',     dest='pkl',           type=str,   action='store', help='extract information from the pickle file (default: None).',                       default=None)
parser.add_argument('--crit-pkl',       metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                          default=None)
parser.add_argument('--bpkl',           metavar='bpkl',    dest='bpkl',          type=str,   action='store', help='enter name of file with all the information about the binodals (default: None).', default=None)
parser.add_argument('--plot-edges',     dest='pe',         action='store_true',  help='plot the edges of the spinodal.')
parser.add_argument('--plot-crits',     dest='pc',         action='store_true',  help='plot the critical points.'      )
parser.add_argument('--no-rtw',         dest='nrtw',       action='store_true',  default=False, help="Dont print out the runtime warning."    )
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

	print("Plotted out second sweep!", flush=True)

	# extract island information
	f = open(args.pkl, 'rb')
	islands = pickle.load(f)
	f.close()

	print(f"Number of islands is {len(islands)}...")

	if len(islands) == 1:
		print(f"Tangent tracing ought to do the job.")
		P.tangent_tracing(ax)

	elif len(islands) == 2:
		pass 

	elif len(islands) == 4:

		# transform the islands and get the hulls
		print(f"Getting the hulls for the islands...", flush=True, end=' ')
		hull_paths = transform_islands(islands) 
		print("done!", flush=True)

		# find the island in the middle
		# the island whose center is closest to [1/3, 1/3] is the central island. 
		print(f"Get the index of the central island...", flush=True, end=' ')
		c_idx = find_central_island(islands)
		print
		print(f"c_idx = {c_idx} and done!", flush=True)

		try:
			print(f"Extracting binodal information from a pickled file...", flush=True, end=' ')
			f = open(args.bpkl, 'rb')
			b_info = pickle.load(f)
			f.close()
			binodals = b_info[0]
			keys_to_center = b_info[1]
			curves_inside_center = b_info[2]
			print("done!", flush=True)

		except:
			print(f"Getting the binodals and other stuff...", flush=True, end=' ')
			f = open(args.bpkl, 'wb')
			binodals              = dict() 
			keys_to_center        = []
			curves_inside_center  = []

			for i in range(len(islands)):
				for j in range(i+1, len(islands)):
					results = P.sym_mu_ps.perform_sweep(islands[i], islands[j])

					# perform the sweep 
					sol1, sol2 = P.sym_mu_ps.binodal_finder(results[0], results[1], hull_paths[i], hull_paths[j])

					# binodals[(i,j)] = [sol1, sol2]
					# if c_idx in (i,j):
					# 	keys_to_center.append((i,j))
					# binodals[(central_island, outside_island)] = [central_binodal, outside_binodal]
					if i == c_idx:
						curves_inside_center.append(sol1)
						binodals[(i,j)] = [sol1, sol2]
						keys_to_center.append((i,j))
					else:
						curves_inside_center.append(sol2)
						binodals[(j, i)] = [sol2, sol1]
						keys_to_center.append((j,i))
					
			
			b_info = [binodals, keys_to_center, curves_inside_center]
			pickle.dump(b_info, f)
			f.close()
			print("done!", flush=True)

		# for key in binodals:
		# 	ax.scatter(binodals[key][0][:,0], 1-binodals[key][0][:,0]-binodals[key][0][:,1], binodals[key][0][:,1], c='gold', s=0.5)
		# 	ax.scatter(binodals[key][1][:,0], 1-binodals[key][1][:,0]-binodals[key][1][:,1], binodals[key][1][:,1], c='pink', s=0.5)

		# start curating curves
		# sort curves using theta 
		# to do so, find the center of the curves

		center = np.array([1/3,1/3,1/3]) 
		print(f"center of central curves is {center}.", flush=True) 
		
		# now, start sorting wrt theta
		linestrings  = []
		for c_idx, c_curve in enumerate(curves_inside_center):
			# c_curve is a 3d array 
			angles = np.arctan2(c_curve[:,1] - center[1], c_curve[:,0] - center[0])

			# sorting curve 
			c_curve = c_curve[np.argsort(angles)]
			binodals[keys_to_center[c_idx]][0] = binodals[keys_to_center[c_idx]][0][np.argsort(angles)]
			binodals[keys_to_center[c_idx]][1] = binodals[keys_to_center[c_idx]][1][np.argsort(angles)]

			# find the biggest angular difference 
			angles = angles[np.argsort(angles)]
			max_delta_theta = np.max(angles[1:] - angles[:-1])

			if max_delta_theta > np.pi/2:
				dists        = np.linalg.norm(c_curve[1:] - c_curve[:-1], axis=1)
				max_dist_idx = np.argmax(dists)
				c_curve      = np.vstack((c_curve[max_dist_idx+1:], c_curve[:max_dist_idx+1]))
				binodals[keys_to_center[c_idx]][0] = np.vstack((binodals[keys_to_center[c_idx]][0][max_dist_idx+1:], binodals[keys_to_center[c_idx]][0][:max_dist_idx+1]))
				binodals[keys_to_center[c_idx]][1] = np.vstack((binodals[keys_to_center[c_idx]][1][max_dist_idx+1:], binodals[keys_to_center[c_idx]][1][:max_dist_idx+1]))

			curves_inside_center[c_idx] = c_curve
			# print(curves_inside_center[c_idx])

			linestrings.append(LineString(curves_inside_center[c_idx][:,0:2]))

		ls_len = len(linestrings)

		# print(f"len(curves_inside_center) = {len(curves_inside_center)}.")
		# cols = ["white", "pink", "gold"]
		# for cidx, curve in enumerate(curves_inside_center):
		# 	ax.plot(curve[:,0], 1-curve[:,0]-curve[:,1], curve[:,1], c=cols[cidx], lw=0.5)
		
		for i in range(ls_len):
			for j in range(i+1, ls_len):
				# print(f"i = {i}, j = {j}...", flush=True)
				# figure out which curve comes first, per theta 
				d1 = (curves_inside_center[i][0][0:2] - center[0:2])/np.linalg.norm(curves_inside_center[i][0][0:2] - center[0:2])
				d2 = (curves_inside_center[j][0][0:2] - center[0:2])/np.linalg.norm(curves_inside_center[j][0][0:2] - center[0:2])
				# print(f"curves_inside_center[i][0:10, 0:2] = {curves_inside_center[i][0:10, 0:2]}, curves_inside_center[j][0][0:2] = {curves_inside_center[j][0:10,0:2]}...")
				# print(f"d1 = {d1}, d2 = {d2}...")
				t =  np.cross(d1[0:2], d2[0:2]) 
				# print(f"t = {t}")
				if t > 0:
					# d2 comes before d1
					intrsxn = linestrings[i].intersection(linestrings[j])
					# print(f"intrsxn = {intrsxn}...")
					try:
						p       = np.array([intrsxn.x, intrsxn.y])
						ax.scatter(p[0], 1-p[0]-p[1], p[1], c='black', s=2)
					except:
						p       = np.array([intrsxn.geoms[0].x, intrsxn.geoms[0].y])
						for point in intrsxn.geoms:
							ax.scatter(point.x, 1-point.y-point.x, point.y, c='black', s=2)
					# kill everything before p for d2, everything after p for d1

					axis = (p-center[0:2])/np.linalg.norm(p-center[0:2])
					direction = (curves_inside_center[i][:,0:2]-center[0:2])/np.linalg.norm(curves_inside_center[i][:,0:2]-center[0:2], axis=1)[:, np.newaxis]

					kept = np.cross(direction[:,0:2], axis)>=0 # direction[np.cross(direction[:,0:2], axis)>=0]
					curves_inside_center[i] = curves_inside_center[i][kept]
					binodals[keys_to_center[i]][0] = binodals[keys_to_center[i]][0][kept]
					binodals[keys_to_center[i]][1] = binodals[keys_to_center[i]][1][kept]

					direction = (curves_inside_center[j][:,0:2]-center[0:2])/np.linalg.norm(curves_inside_center[j][:,0:2]-center[0:2], axis=1)[:, np.newaxis]
					kept = np.cross(direction[:,0:2], axis)<0 # direction[np.cross(direction[:,0:2], axis)<0]
					curves_inside_center[j] = curves_inside_center[j][kept]
					binodals[keys_to_center[j]][0] = binodals[keys_to_center[j]][0][kept]
					binodals[keys_to_center[j]][1] = binodals[keys_to_center[j]][1][kept]					
					# print(f"curves_inside_center[i][0:10, 0:2] = {curves_inside_center[i][0:10, 0:2]}, \ncurves_inside_center[j][0][0:2] = {curves_inside_center[j][0:10,0:2]}...")
				
				else:
					# d2 comes before d1
					intrsxn = linestrings[i].intersection(linestrings[j])
					# print(f"intrsxn = {intrsxn}...")
					try:
						p       = np.array([intrsxn.x, intrsxn.y])
						ax.scatter(p[0], 1-p[0]-p[1], p[1], c='black', s=2)
					except:
						p       = np.array([intrsxn.geoms[0].x, intrsxn.geoms[0].y])
						for point in intrsxn.geoms:
							ax.scatter(point.x, 1-point.x-point.y, point.y, c='black', s=2)
					# kill everything before p for d2, everything after p for d1
					axis = (p-center[0:2])/np.linalg.norm(p-center[0:2])
					direction = (curves_inside_center[j][:,0:2]-center[0:2])/np.linalg.norm(curves_inside_center[j][:,0:2]-center[0:2], axis=1)[:, np.newaxis]
					kept = np.cross(direction[:,0:2], axis)>=0 # direction[np.cross(direction[:,0:2], axis)>=0]
					curves_inside_center[j] = curves_inside_center[j][kept]
					binodals[keys_to_center[j]][0] = binodals[keys_to_center[j]][0][kept]
					binodals[keys_to_center[j]][1] = binodals[keys_to_center[j]][1][kept]	

					direction = (curves_inside_center[i][:,0:2]-center[0:2])/np.linalg.norm(curves_inside_center[i][:,0:2]-center[0:2], axis=1)[:, np.newaxis]
					kept = np.cross(direction[:,0:2], axis[0:2])<0 # direction[np.cross(direction[:,0:2], axis[0:2])<0]
					curves_inside_center[i] = curves_inside_center[i][kept]
					binodals[keys_to_center[i]][0] = binodals[keys_to_center[i]][0][kept]
					binodals[keys_to_center[i]][1] = binodals[keys_to_center[i]][1][kept]
		
		# go to an island, and iterrogate all the binodals there. 
		# whichever binodals intersect, use them to construct the tie points
		bkeys = list(binodals.keys())

		for i in range(len(bkeys)):
			for j in range(i+1, len(bkeys)):
				# there is an island in which the two sets of binodals have an arm 
				# see if they intersect. to do this, convert binodals to linestrings	
				if bkeys[i][0] in bkeys[j]:
					L1 = LineString(binodals[bkeys[i]][0][:,0:2])
					print(f"L1 = {binodals[bkeys[i]][0][:10,0:2]}")
					L2 = LineString(binodals[bkeys[j]][bkeys[j].index(bkeys[i][0])][:,0:2])
					print(f"L2 = {binodals[bkeys[j]][bkeys[j].index(bkeys[i][0])][:10,0:2]}")
					intr = L1.intersection(L2)
					print(f"For i = {bkeys[i]}, j = {bkeys[j]}, intersection = {intr}.", flush=True)
					p       = np.array([intr.x, intr.y])
					ax.scatter(intrsxn.x, 1-intrsxn.x-intrsxn.y, intrsxn.y, c='peru', s=2)

				elif bkeys[i][1] in bkeys[j]:
					L1 = LineString(binodals[bkeys[i]][1][:,0:2])
					L2 = LineString(binodals[bkeys[j]][bkeys[j].index(bkeys[i][1])][:,0:2])
					intr = L1.intersection(L2)
					print(f"For i = {bkeys[i]}, j = {bkeys[j]}, intersection = {intr}.", flush=True)
					
					ax.scatter(intr.x, 1-intr.x-intr.y, intr.y, c='peru', s=2)

				else:
					print(f"No shared islands for i = {bkeys[i]}, j = {bkeys[j]}.", flush=True)					
					# fuck em
					pass 


		# print(f"len(curves_inside_center) = {len(curves_inside_center)}.")
		cols = ["white", "pink", "gold", "grey", "skyblue", "lavender"]
		idx = 0
		for bkey in binodals:
			ax.scatter(binodals[bkey][0][:,0], 1-binodals[bkey][0][:,0]-binodals[bkey][0][:,1], binodals[bkey][0][:,1], s=0.5, c=cols[idx%len(cols)])
			ax.scatter(binodals[bkey][1][:,0], 1-binodals[bkey][1][:,0]-binodals[bkey][1][:,1], binodals[bkey][1][:,1], s=0.5, c=cols[idx%len(cols)])
			idx += 1
		# for cidx, curve in enumerate(curves_inside_center):
		# 	ax.scatter(curve[:,0], 1-curve[:,0]-curve[:,1], curve[:,1], c=cols[cidx], s=0.5)

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

