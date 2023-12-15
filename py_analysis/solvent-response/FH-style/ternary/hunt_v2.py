import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from matplotlib.path import Path
from scipy.optimize import fsolve
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
from scipy.spatial.distance import cdist
import sys
import argparse
import linecache
import mpltern
np.set_printoptions(threshold=sys.maxsize)
import warnings 
import tangent
import ternary
import time
from sklearn.cluster import DBSCAN
import pickle 
import copy
import argparse 
from scipy.spatial import ConvexHull

parser = argparse.ArgumentParser(description="Locate the /c/ritical points on the /spin/odal diagram. This program will create one plot, and you can customize what you want on the plot.")
parser.add_argument('--chisc',  metavar='chi_sc',  dest='chi_sc',  type=float,   action='store', help='enter A-C exchange parameter.' )
parser.add_argument('--chips',  metavar='chi_ps',  dest='chi_ps',  type=float,   action='store', help='enter A-B exchange parameter.' )
parser.add_argument('--chipc',  metavar='chi_pc',  dest='chi_pc',  type=float,   action='store', help='enter B-C exchange parameter.' )
parser.add_argument('-vs',      metavar='vs',      dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',      dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',      dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--search-density',      dest='sd',        type=int,         action='store', help='search density of hunt.')
parser.add_argument('--ternary',             dest='ternary',   action='store_true',  default=False, help='make the output a ternary plot.')
parser.add_argument('--no-rtw',              dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--pkl',                 dest='pkl',       action='store', type=str,  help='extract information from the pickle file.')
parser.add_argument('--img-name',            dest='img',       action='store', type=str,  default="None", help='name of the image to be created (default: all of the inputs in the imagename).')
args = parser.parse_args()

def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f"beep.\n"
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

def perform_sweep(island1, island2):

	print(f"Performing a sweep to check for chemical potentials...", flush=True)
	chem_pot_a_upper = mu_a(island1[:,0], island1[:,1])
	chem_pot_b_upper = mu_b(island1[:,0], island1[:,1])
	chem_pot_c_upper = mu_c(island1[:,0], island1[:,1])

	masks = np.isinf(chem_pot_a_upper) | np.isnan(chem_pot_a_upper) | np.isinf(chem_pot_b_upper) |\
	 np.isnan(chem_pot_b_upper) | np.isinf (chem_pot_c_upper) | np.isnan (chem_pot_c_upper)
	
	chem_pot_a_upper = chem_pot_a_upper [~masks]
	chem_pot_b_upper = chem_pot_b_upper [~masks]
	chem_pot_c_upper = chem_pot_c_upper [~masks]

	# get only the relevant phi_uppers
	island1 = island1[~masks]
	# print(island1)

	chem_pot_a_lower = mu_a(island2[:,0], island2[:,1])
	chem_pot_b_lower = mu_b(island2[:,0], island2[:,1])
	chem_pot_c_lower = mu_c(island2[:,0], island2[:,1])

	# create a mask for all the good lower chemical potentials
	masks = np.isinf(chem_pot_a_lower) | np.isnan(chem_pot_a_lower) | np.isinf(chem_pot_b_lower) |\
	np.isnan(chem_pot_b_lower) | np.isinf (chem_pot_c_lower) | np.isnan (chem_pot_c_lower)

	# get only the relevant phi_lowers
	chem_pot_a_lower = chem_pot_a_lower [~masks]
	chem_pot_b_lower = chem_pot_b_lower [~masks]
	chem_pot_c_lower = chem_pot_c_lower [~masks]

	# get only the relevant phi_lowers
	island2 = island2[~masks]

	# get an array of these values
	mu_upper   = np.array ([chem_pot_a_upper, chem_pot_b_upper, chem_pot_c_upper]).T
	mu_lower   = np.array ([chem_pot_a_lower, chem_pot_b_lower, chem_pot_c_lower]).T

	# calculate all the distances
	distances       = np.linalg.norm (mu_upper[:, np.newaxis] - mu_lower, axis=2)
	print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
	closest_indices = np.argmin(distances, axis=1)

	# partition, partition... 
	island2         = island2[closest_indices]
	mu_lower        = mu_lower[closest_indices]
	min_distances   = distances[np.arange(len(mu_upper)), closest_indices]

	# island2 has a bunch of repeated rows 
	# closest_indices has a bunch of repeated indices. The index at which the repeated indices are present
	# give us the indices of the location of the thing in A.
	# find the unique indices in closest_indices
	# print(f"closest indices = {closest_indices[0:10]}", flush=True)
	# print(f"island1 = {island1}")
	'''
	unique_indices = np.unique(closest_indices)
	# print(f"len(unique_indices) = {len(unique_indices)}")
	
	res = []
	for index in unique_indices:
		# find rows in island1 which are closest to this particular unique index
		rows_matching_index = np.where(closest_indices == index)[0]
		# print(f"rows_matching_index = {rows_matching_index}")
		# find the gaps in chemical potential between the points in island1 and the one point in island2
		MU_A = mu_a(island1[rows_matching_index][:,0], island1[rows_matching_index][:,1]) - mu_a(island2[index][0], island2[index][1])
		MU_B = mu_b(island1[rows_matching_index][:,0], island1[rows_matching_index][:,1]) - mu_b(island2[index][0], island2[index][1])
		MU_C = mu_c(island1[rows_matching_index][:,0], island1[rows_matching_index][:,1]) - mu_c(island2[index][0], island2[index][1])
		DMU  = np.array([MU_A, MU_B, MU_C]).T
		distances = np.linalg.norm(DMU, axis=1)

		min_distance_index = rows_matching_index[np.argmin(distances)]
		res.append(island1[min_distance_index])
	
	island1 = np.array(res)
	island2 = island2[unique_indices]
	''' 
	return [island1, island2, min_distances]


def binodal_finder(island1, island2, hull1, hull2, vs, vc, vp, chi_ps, chi_pc, chi_sc):

	sol1 = np.empty((0,3))
	sol2 = np.empty((0,3))

	for idx, i1 in enumerate(island1):
		def mu_equations(phi):
			eq1 = mu_a(phi[0], island1[idx,1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], island1[idx,1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], island1[idx,1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root = fsolve(mu_equations, [island1[idx,0], island2[idx,0], island2[idx,1]])
		p1   = np.array([root[0], island1[idx,1], 1-root[0]-island1[idx,1]])
		p2   = np.array([root[1], root[2], 1-root[1]-root[2]])

		if (np.abs(np.array(mu_equations(root)))>1e-6).any():
			continue
		else:
			if (ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0) or (ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0):
				continue 

			elif np.isnan(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				continue 
			
			elif np.isinf(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				continue 

			elif np.linalg.norm(p1-p2) < 1e-4:
				continue

			#elif (not hull1.contains_point(p1[0:2])) or (not hull2.contains_point(p2[0:2])):
			#	continue 

			else:
				sol1 = np.vstack((sol1, p1))
				sol2 = np.vstack((sol2, p2))


		def mu_equations(phi):
			eq1 = mu_a(island1[idx,0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(island1[idx,0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(island1[idx,0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		root = fsolve(mu_equations, [island1[idx,1], island2[idx,0], island2[idx,1]])
		p1   = np.array([island1[idx,0], root[0], 1-island1[idx,0]-root[0]])
		p2   = np.array([root[1], root[2], 1-root[1]-root[2]])

		if (np.abs(np.array(mu_equations(root)))>1e-6).any():
			continue
		else:
			if (ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0) or (ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0):
				continue 

			elif np.isnan(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				continue 
			
			elif np.isinf(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
				continue

			elif np.linalg.norm(p1-p2) < 1e-4:
				continue

			#elif (not hull1.contains_point(p1[0:2])) or (not hull2.contains_point(p2[0:2])):
			#	continue 

			else:
				sol1 = np.vstack((sol1, p1))
				sol2 = np.vstack((sol2, p2))
		

	return sol1, sol2 

if __name__=="__main__":

	print ("Revving up the program...", flush=True)

	####################################################
	# set up the inputs
	lsize = 3
	font = {'color':  'black',
		'weight': 'normal',
		'size': lsize}

	fig = plt.figure(figsize=(8,8))
	if args.ternary:
		ax = fig.add_subplot(projection="ternary")
	else:
		ax = plt.axes()

	chi_sc = args.chi_sc
	chi_ps = args.chi_ps
	chi_pc = args.chi_pc
	vs     = args.vs
	vc     = args.vc
	vp     = args.vp
	####################################################

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

	######################################################
	# FIND PHI_B GIVEN PHI_A
	mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - vs/vp * phi_b - vs/vc * (1-phi_a-phi_b) + vs * (phi_b**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_sc + phi_b * (1-phi_a-phi_b) * (chi_ps + chi_sc - chi_pc) ) 
	mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vp/vs * phi_a - vp/vc * (1-phi_a-phi_b) + vp * (phi_a**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_pc + phi_a * (1-phi_a-phi_b) * (chi_ps + chi_pc - chi_sc) )
	mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/vs * phi_a - vc/vp * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )
	######################################################

	# get the island file
	f = open(args.pkl, 'rb')

	# extract information from that file
	islands = pickle.load(f)

	# close the file
	f.close()

	# plot the islands out 
	cols  = ["coral", "steelblue", "forestgreen", "lavender"]

	# pad the islands
	for idx, island in enumerate(islands):
		islands[idx] = np.array(np.vstack([0.001+(0.999-0.001)/(args.sd-1)*island[:,1], 0.001+(0.999-0.001)/(args.sd-1)*island[:,0]])).T 
		xllim = np.min(islands[idx][:,0])*0.95
		xulim = np.max(islands[idx][:,0])*1.05 if np.max(islands[idx][:,0])*1.05 < 1 else 1
		yllim = np.min(islands[idx][:,1])*0.95
		yulim = np.max(islands[idx][:,1])*1.05 if np.max(islands[idx][:,1]*1.05) < 1 else 1
		print(f"islands[{idx}].shape = {islands[idx].shape}")

		if idx == 2:					
			x = np.linspace(xllim, xulim, 500)
			y = np.linspace(yllim, yulim, 500)
			xv, yv = np.meshgrid(x, y)
			to_keep = (ternary.stab_crit(xv, yv, vs, vc, vp, chi_ps, chi_pc, chi_sc))>=0
			x_keep = xv[to_keep]
			y_keep = yv[to_keep]
			islands[idx] = np.vstack((islands[idx], np.array([x_keep, y_keep]).T)) 
			print(f"islands[{idx}].shape = {islands[idx].shape}")

		print(f"islands[{idx}] = {islands[idx][0:10]}", flush=True)
	exit()
	hull_paths = []
	for idx, island in enumerate(islands):
		# print(f"Island = \n{island}")
		ax.scatter(island[:,0], island[:,1], c=cols[idx%4], s=1)
		
		hull = ConvexHull(islands[idx])
		for simplex in hull.simplices:
			# print(f"hull = {islands[idx][simplex]}")
			ax.plot(islands[idx][simplex,0], islands[idx][simplex,1], 'k-')
		hull_paths.append( Path(islands[idx][hull.vertices]) )
	print("Done transforming the islands.", flush=True)

	# run a island by island sweep
	print(f"len(islands) = {len(islands)}")
	L = len(islands)
	for i in [0, 1, 3]: # range(L):
		for j in [2]: # range(i+1, L):
			# obtain nearest chem_pots
			print(f"i={i}, j={j}", flush=True)
			results = perform_sweep(islands[i], islands[j])
			print(f"results[0] = \n{results[0][:10]}", flush=True)
			print(f"results[1] = \n{results[1][:10]}", flush=True)
			# perform the sweep
			sol1, sol2 = binodal_finder(results[0], results[1], hull_paths[i], hull_paths[j], vs, vc, vp, chi_ps, chi_pc, chi_sc)
			ax.scatter(sol1[:,0], sol1[:,1], c="black",  s=0.5)
			ax.scatter(sol2[:,0], sol2[:,1], c="crimson", s=0.5)
			print(f"sol1 = \n{sol1[:10]}")
			print(f"sol2 = \n{sol2[:10]}")

			for k in range(0, len(sol2), 100):
				ax.plot([sol1[k,0],sol2[k,0]], [sol1[k,1], sol2[k,1]], lw=0.25, ls='--', c='pink', alpha=0.5)

	print(f"len(islands) = {len(islands)}")

	for i in [1]: # range(L):
		for j in [0, 2, 3]: # range(i+1,L): # [0,2,3]:
			# obtain nearest chem_pots
			print(f"i={L-(i+1)}, j={L-(j+1)}", flush=True)
			results = perform_sweep(islands[L-(i+1)], islands[L-(j+1)])
			print(f"results[0] = \n{results[0][:10]}", flush=True)
			print(f"results[1] = \n{results[1][:10]}", flush=True)
			# perform the sweep
			sol1, sol2 = binodal_finder(results[0], results[1], hull_paths[L-(i+1)], hull_paths[L-(j+1)], vs, vc, vp, chi_ps, chi_pc, chi_sc)
			ax.scatter(sol1[:,0], sol1[:,1], c="black",  s=0.5)
			ax.scatter(sol2[:,0], sol2[:,1], c="crimson", s=0.5)
			print(f"sol1 = \n{sol1[:10]}")
			print(f"sol2 = \n{sol2[:10]}")

			for k in range(0, len(sol2)):
				ax.plot([sol1[k,0],sol2[k,0]], [sol1[k,1], sol2[k,1]], lw=0.25, ls='--', c='pink', alpha=0.5)


	# results = perform_sweep(islands[0], islands[2])

	# ax.scatter(results[0][:,0], results[0][:,1], c='gray')
	# ax.scatter(results[1][:,0], results[1][:,1], c='orange')

	# perform the sweep
	# sol1, sol2 = binodal_finder(results[0], results[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)

	# ax.scatter(sol1[:,0], sol1[:,1], c="black",  s=0.5)
	# ax.scatter(sol2[:,0], sol2[:,1], c="crimson", s=0.5)

	# for i in range(len(sol2)):
	# 	ax.plot([sol1[i,0],sol2[i,0]], [sol1[i,1], sol2[i,1]], lw=0.25, ls='--', c='pink')


	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		plt.savefig (img_name, dpi=1200, bbox_inches="tight")
	elif args.ternary:
		plt.savefig (f"signs_tern-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)
	else:
		plt.savefig (f"signs_reg-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)

	print ("Completed heat map computation.")

