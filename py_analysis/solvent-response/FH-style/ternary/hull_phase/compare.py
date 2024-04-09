import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import FreeEnergy as FE
import linecache 
import pickle
import warnings 
import time
import mu
import random
import mpltern
import argparse
import ternary
from scipy.optimize import fsolve
from itertools import combinations
from scipy.spatial import ConvexHull

parser = argparse.ArgumentParser(description='Check the triplets that form a good solution.')
parser.add_argument('--hull' , metavar='hull',  dest='hull',  type=str, action='store', help='location of serialized hull structure.')
parser.add_argument('--tdump', metavar='tdump', dest='tdump', type=str, action='store', help='location of depositing good triplets.' )
parser.add_argument('--tget',  metavar='tget',  dest='tget',  type=str, action='store', help='location of obtaining triplets.'       )
args = parser.parse_args()

###################################################
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
###################################################

def min_sum_mutual_distances(points):
    """
    Find three points from an array of points such that the sum of mutual distances is minimum.

    Parameters:
        points (numpy.ndarray): Array of shape (N, 3) containing coordinates of points.

    Returns:
        tuple: Indices of the three points.
    """
    min_sum = float('inf')
    min_indices = None

    # Generate combinations of three points
    for combo in combinations(range(len(points)), 3):
        # Calculate sum of mutual distances
        sum_distances = sum(np.linalg.norm(points[combo[i]] - points[combo[j]]) for i in range(3) for j in range(i+1, 3))
        
        # Update minimum sum and indices if current sum is smaller
        if sum_distances < min_sum:
            min_sum = sum_distances
            min_indices = combo

    return min_indices

def custom_sort(arrays):
    # Convert the list of arrays into a NumPy array
    arrays = np.array(arrays)
    # Sort the arrays based on the first element
    sorted_indices = np.argsort(arrays[:, 0]) 
    sorted_arrays = arrays[sorted_indices]

    # Check for elements with first elements within 1e-6 of each other
    first_elements = sorted_arrays[:, 0]
    diff = np.abs(np.diff(first_elements))
    split_indices = np.where(diff > 1e-6)[0] + 1 

    # Sort the subarrays based on the second element
    final_sorted_arrays = np.split(sorted_arrays, split_indices)
    final_sorted_arrays = [np.array(sorted(subarr, key=lambda x: x[1])) for subarr in final_sorted_arrays]

    # Concatenate the sorted subarrays
    sorted_arrays = np.concatenate(final_sorted_arrays)

    return sorted_arrays[:,0:2]

def triangle_area(p1, p2, p3):
	area = 0.5 * np.abs(np.cross(p2-p1, p3-p1))
	return area 

if __name__=="__main__":

	start = time.time()

	hull_file = open(args.hull, 'rb')
	F = pickle.load(hull_file)
	hull_file.close() 

	inputs = dict()
	inputs["vs"] = F.vs
	inputs["vc"] = F.vc
	inputs["vp"] = F.vp 
	inputs["chi_sc"] = F.chi_sc
	inputs["chi_ps"] = F.chi_ps 
	inputs["chi_pc"] = F.chi_pc

	df   = pd.read_csv(args.tget, sep='\|', engine='python', names=["vs", "vc", "vp", "chisc", "chips", "chipc", "phi_s", "phi_p", "l0", "l1", "l2", "phi_s1", "phi_p1", "phi_s2", "phi_p2", "phi_s3", "phi_p3", "w1", "w2", "w3"], skiprows=1)
	phi_s1, phi_s2, phi_s3 = df["phi_s1"].values, df["phi_s2"].values, df["phi_s3"].values
	phi_p1, phi_p2, phi_p3 = df["phi_p1"].values, df["phi_p2"].values, df["phi_p3"].values

	fig = plt.figure(figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	ax.set_tlabel ("$\\phi _{S}$") # ('Vol. frac. A')
	ax.set_llabel ("$\\phi _{C}$") # ('Vol. frac. C')
	ax.set_rlabel ("$\\phi _{P}$") # ('Vol. frac. B')
	ax.set_tlim(0, 1)
	ax.set_llim(0, 1)
	ax.set_rlim(0, 1)
	positions = ['tick1', 'tick2']
	for position in positions:
		ax.taxis.set_ticks_position(position)
		ax.laxis.set_ticks_position(position)
		ax.raxis.set_ticks_position(position)

	ax.taxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.raxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.laxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])

	ax.scatter(F.crits[:,0], 1-F.crits[:,0]-F.crits[:,1], F.crits[:,1], c='black', s=0.9)
	M    = mu.sym_mu_ps(inputs, F.spinodal)

	def delta_mu(P):
		eq1 = M.delta_mu_s(P[0], P[1], P[4], P[5])
		eq2 = M.delta_mu_p(P[0], P[1], P[4], P[5])
		eq3 = M.delta_mu_c(P[0], P[1], P[4], P[5])
		eq4 = M.delta_mu_s(P[2], P[3], P[4], P[5])
		eq5 = M.delta_mu_p(P[2], P[3], P[4], P[5])
		eq6 = M.delta_mu_c(P[2], P[3], P[4], P[5])
		return [eq1, eq2, eq3, eq4, eq5, eq6]

	dump = open(args.tdump, 'w')
	for i in range(len(phi_s1)):
		phi_tests = [phi_s1[i], phi_p1[i], phi_s2[i], phi_p2[i], phi_s3[i], phi_p3[i]]
		root = fsolve(delta_mu, phi_tests)
		p1   = np.array([root[0], root[1]])
		p2   = np.array([root[2], root[3]])
		p3   = np.array([root[4], root[5]])
		area = triangle_area(p1, p2, p3)
		if (np.abs(delta_mu(root)) > 1e-6).any():
			continue
		elif np.linalg.norm(p1-p2)<1e-6 or np.linalg.norm(p1-p3)<1e-6 or np.linalg.norm(p2-p3)<1e-6 or area < 1e-3:
			continue
		else:
			p = np.array([p1, p2, p3])
			p = custom_sort(p)
			print(f"Good root = {root}", flush=True)
			print(f"p1 = {p1}, p2 = {p2}, p3 = {p3} | error = {delta_mu(root)}")
			print(f"area = {area}")
			ax.plot([p[0][0],p[1][0],p[2][0],p[0][0]], [1-p[0][0]-p[0][1], 1-p[1][0]-p[1][1], 1-p[2][0]-p[2][1], 1-p[0][0]-p[0][1]], [p[0][1], p[1][1], p[2][1], p[0][1]], lw=0.25, c='steelblue')
			wts = np.random.uniform(low=0, high=1, size=3)
			wts = wts/np.sum(wts)
			pt  = wts[0]*p[0] + wts[1]*p[1] + wts[2]*p[2]			
			dump.write(f"{F.vs} | {F.vc} | {F.vp} | {F.chi_sc} | {F.chi_ps} | {F.chi_pc} | {pt[0]} | {pt[1]} | 0 | 0 | 1 | {p[0][0]} | {p[0][1]} | {p[1][0]} | {p[1][1]} | {p[2][0]} | {p[2][1]} | {wts[0]} | {wts[1]} | {wts[2]}\n")
			
	dump.close()


	df   = pd.read_csv("double.db", sep='\|', engine='python', names=["vs", "vc", "vp", "chisc", "chips", "chipc", "phi_s", "phi_p", "l0", "l1", "l2", "phi_s1", "phi_p1", "phi_s2", "phi_p2", "phi_s3", "phi_p3", "w1", "w2", "w3"], skiprows=1)
	phi_s1, phi_s2 = df["phi_s1"].values, df["phi_s2"].values
	phi_p1, phi_p2 = df["phi_p1"].values, df["phi_p2"].values
	def delta_mu(P, phi_s):
		eq1 = M.delta_mu_s(phi_s, P[0], P[1], P[2])
		eq2 = M.delta_mu_p(phi_s, P[0], P[1], P[2])
		eq3 = M.delta_mu_c(phi_s, P[0], P[1], P[2])
		return [eq1, eq2, eq3]

	print(f"Entering two phase calcs...")
	for i in range(len(phi_s1)):
		phi_tests = [phi_p1[i], phi_s2[i], phi_p2[i]]
		root = fsolve(delta_mu, phi_tests, args=(phi_s1[i]))
		p1   = np.array([phi_s1[i], root[0]])
		p2   = np.array([root[1], root[2]])
		if (np.abs(delta_mu(root, phi_s1[i])) > 1e-6).any():
			continue
		elif np.linalg.norm(p1-p2)<1e-2:
			continue
		else:
			p = np.array([p1, p2])
			p = custom_sort(p)
			print(f"Good root = {root}", flush=True)
			print(f"p1 = {p1}, p2 = {p2}| error = {delta_mu(root, phi_s1[i])}")
			ax.plot([p[0][0],p[1][0]], [1-p[0][0]-p[0][1], 1-p[1][0]-p[1][1]], [p[0][1], p[1][1]], c='coral', lw=0.5, ls='--', alpha=0.1)

	fig.savefig("net.png", dpi=1000, bbox_inches="tight")

	stop = time.time()
	print(f"Time for database creation is {stop-start} seconds.", flush=True)
