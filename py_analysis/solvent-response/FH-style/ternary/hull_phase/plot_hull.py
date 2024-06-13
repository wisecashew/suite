import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpltern
import ternary
from scipy.optimize import fsolve
import FreeEnergy as FE
import pickle
import argparse 
import mu

parser = argparse.ArgumentParser(description='')
parser.add_argument('--hull', metavar='hull', dest='hull', type=str, action='store', help='location of serialized hull structure.')
parser.add_argument('--img',  metavar='img',  dest='img',  type=str, action='store', help='name of image to be created (default: my_database).', default="my_database")
args = parser.parse_args()

def custom_sort(arrays):
	# Convert the list of arrays into a NumPy array
	arrays = np.array(arrays)

	# Sort the arrays based on the first element
	sorted_indices = np.argsort(arrays[:, 0]) 
	sorted_arrays = arrays[sorted_indices]

	# Check for elements with first elements within 1e-6 of each other
	first_elements = sorted_arrays[:, 0]
	diff           = np.abs(np.diff(first_elements))
	split_indices  = np.where(diff > 1e-6)[0] + 1 

	# Sort the subarrays based on the second element
	final_sorted_arrays = np.split(sorted_arrays, split_indices)
	final_sorted_arrays = [np.array(sorted(subarr, key=lambda x: x[1])) for subarr in final_sorted_arrays]

	# Concatenate the sorted subarrays
	sorted_arrays = np.concatenate(final_sorted_arrays)

	return sorted_arrays[:,0:2]

if __name__=="__main__":

	fig = plt.figure(figsize=(10,10))
	ax  = fig.add_subplot(projection="ternary")

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

	nsimplices = len(F.simplices)
	M = mu.sym_mu_ps(inputs, F.spinodal)
	one_phase_idx   = np.arange(nsimplices)[F.num_phases==1]
	two_phase_idx   = np.arange(nsimplices)[F.num_phases==2]
	three_phase_idx = np.arange(nsimplices)[F.num_phases==3]

	two_phase_idx   = two_phase_idx[::100]
	for tdx, tpi in enumerate(two_phase_idx):
		print(f"@ {tdx}/{len(two_phase_idx)}...")
		simplex      = F.surface[F.simplices[tpi]]
		vertices     = F.from_xy_to_phi(simplex)
		weights      = np.random.uniform(low=0, high=1, size=3)
		weights      = weights/np.sum(weights) 
		point        = np.sum(weights.reshape(-1,1)*vertices, axis=0)
		results      = F.calc_phase_split_specific(point, tpi)
		phase_comps  = custom_sort(results[1][0:2, 0:2])
		p1, p2       = phase_comps[0], phase_comps[1]
		ax.scatter(p1[0], 1-p1[0]-p1[1], p1[1], c='coral', s=0.5)
		ax.scatter(p2[0], 1-p2[0]-p2[1], p2[1], c='steelblue', s=0.5)

	fig.savefig(args.img, dpi=1000)
