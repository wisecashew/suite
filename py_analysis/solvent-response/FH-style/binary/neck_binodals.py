import numpy as np 
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root
import sympy as sym
import argparse 
import warnings
import linecache
import time
import cma

EPS    = 1e-4
parser = argparse.ArgumentParser (description="Plots phase diagrams.")
parser.add_argument ("--label",  dest='label', type=str,   action='store', help="Prove the label for the diagram.")
parser.add_argument ("--xlim",  dest='xlim',   type=float, nargs='+', action='store', help="Provide a limits for x-axis.", default=[0,1])
parser.add_argument ("--ylim",  dest='ylim',   type=float, nargs='+', action='store', help="Provide a limits for y-axis.", default=[0,1])
parser.add_argument ("--T",      dest='T',     type=float, nargs='+', action='store', help="Provide a temperature range to plot thing in.", default=None)
parser.add_argument ("--T-loop", dest='Tloop', type=float, action='store', help="Provide a temperature to bisect a loop.", default=None)
parser.add_argument ("--T-UCST", dest='Tucst', type=float, action='store', help="Provide a temperature above which to hunt for a binodal.", default=None)
parser.add_argument ("--T-neck", dest='Tneck', type=float, action='store', help="Provide a temperature above and below which we will hunt for a binodal.", default=None)
parser.add_argument ("--custom-weights", dest='cw',   action='store_true', help="Enter option to draw custom weights.", default=False)
parser.add_argument ("--skip",           dest='skip', type=int, action='store', help="Provide length of jump while looping through points for the binodal.", default=1000)
parser.add_argument ("--vm",     dest='vm',     type=float, nargs='+', action='store', help="Provide vm.",   default=-1)
parser.add_argument ("--vs",     dest='vs',     type=float, nargs='+', action='store', help="Provide vs.",   default=-1)
parser.add_argument ("--xfrmpv", dest='pv',     type=float, nargs='+', action='store', help="Provide PV.",   default=-1)
parser.add_argument ("--xfrmpw", dest='pw',     type=float, nargs='+', action='store', help="Provide PW.",   default=-1)
parser.add_argument ("--emm",    dest="emm",    type=float, nargs='+', action='store', help="Provide EMMA.", default=-1)
parser.add_argument ("--emsa",   dest="emsa",   type=float, nargs='+', action='store', help="Provide EMSA.", default=-1)
parser.add_argument ("--emsn",   dest="emsn",   type=float, nargs='+', action='store', help="Provide EMSN.", default=0)
parser.add_argument ("--ess",    dest="ess",    type=float, nargs='+', action='store', help="Provide ESSA.", default=0)
parser.add_argument ("--csv",    dest='csv',    type=str,   nargs='+', action='store', help='Name of csv file.')
parser.add_argument ("--colors", dest='c',      type=str,   nargs='+', action='store', help='Name of colors.')
parser.add_argument ("--marker", dest='marker', type=str,   nargs='+', action='store', help='String of markers.')
parser.add_argument ("--img",    dest="img",    type=str,   action='store', help="Name of image.", default="bintest.png")
args = parser.parse_args() 

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return ""# f"Probably a math sqrt error (line {lineno}).\n"

warnings.formatwarning = custom_warning_format

def remove_close_rows(array, threshold=1e-6):
	kept_indices   = []
	filtered_array = np.empty ((0,array.shape[0]))
	for i, elem in enumerate(array):
		if i == 0:
			filtered_array = np.vstack((filtered_array, elem))
			kept_indices.append(i)
			continue
		else:
			sieve = (np.linalg.norm(filtered_array - elem, axis=1) < threshold).any()
			if sieve:
				continue
			else:
				filtered_array = np.vstack((filtered_array, elem))
				kept_indices.append(i)

	return filtered_array, np.array(kept_indices)

def delete_close_elements(array):
	"""
	Delete elements in a numpy array whose distance from one another is less than 1e-6,
	but keep the one with the smaller index. Return the pruned array and the indices
	of the kept elements.

	Parameters:
	- array: Numpy array containing elements.

	Returns:
	- pruned_array: Numpy array with close elements removed.
	- kept_indices: Indices of the kept elements.
	"""
	pruned_array = []
	kept_indices = []

	for i in range(len(array)):
		keep = True
		for j in range(i + 1, len(array)):
			# Compute the distance between two elements
			distance = np.linalg.norm(array[i] - array[j])

			# If the distance is less than 1e-6, remove the element with the larger index
			if distance < 1e-6:
				keep = False
				break

		# Keep the element if it's not too close to any other element
		if keep:
			pruned_array.append(array[i])
			kept_indices.append(i)

	return np.array(pruned_array), np.array(kept_indices)

def find_closest_indices_vectorized(A, B):
	"""
	Find the indices of the elements in arrays A and B that have the least distance between them (vectorized).

	Parameters:
	- A: Numpy array of length N
	- B: Numpy array of length M

	Returns:
	- index_A: Index of the closest element in array A
	- index_B: Index of the closest element in array B
	- min_distance: Minimum distance between the closest elements
	"""
	# Compute pairwise distances using broadcasting
	distances = np.abs(A[:, np.newaxis] - B)

	# Find the index of the minimum distance
	index_A, index_B = np.unravel_index(np.argmin(distances), distances.shape)

	# Minimum distance
	min_distance = distances[index_A, index_B]

	return index_A, index_B, min_distance

def index_max_distance(array):
	"""
	Find the index of the element that is maximally distant from the next element in the array.

	Parameters:
	- array: Numpy array

	Returns:
	- index: Index of the element with maximum distance from the next element
	"""
	# Compute absolute differences between consecutive elements
	differences = np.abs(np.diff(array))

	# Find the index of the maximum difference
	index = np.argmax(differences)

	return index

zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa) + (1-pw)*np.exp (-1/T * emsn)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa)/zms(emsa, emsn, pw, T)

def loss_weights(experimental):
	weights = []
	if args.label == "LOOP":
		mask    = experimental[:,1] >= args.Tloop
		phi_top = experimental[:,0][mask]

		mask       = experimental[:,1] <= args.Tloop
		phi_bottom = experimental[:,0][mask]

		for btups in experimental:
			# print(f"btups = {btups}")
			if btups[1] >= args.Tloop:
				phi_diffs = btups[0] - phi_top 
				if len(phi_diffs[phi_diffs<0]) == 0:
					dphi_pos = (1-btups[0])*10
				else:
					dphi_pos = np.min(np.abs(phi_diffs[phi_diffs<0]))

				if len(phi_diffs[phi_diffs>0]) == 0:
					dphi_neg = btups[0]*10 
				else:
					dphi_neg   = np.min(np.abs(phi_diffs[phi_diffs>0]))

			else:
				phi_diffs = btups[0] - phi_bottom 
				if len(phi_diffs[phi_diffs<0]) == 0:
					dphi_pos   = (1-btups[0])*10
				else:
					dphi_pos   = np.min(np.abs(phi_diffs[phi_diffs<0]))

				if len(phi_diffs[phi_diffs>0]) == 0:
					dphi_neg = btups[0]*10
				else:
					dphi_neg = np.min(np.abs(phi_diffs[phi_diffs>0]))
			weights.append(dphi_pos+dphi_neg)
		
		weights = np.array(weights)
		mask = experimental[:,0] <= 0.15
		weights[mask] = 0

	elif args.label == "NECK":
		mask    = experimental[:,1] >= args.Tneck
		phi_top = experimental[:,0][mask]

		mask       = experimental[:,1] <= args.Tneck
		phi_bottom = experimental[:,0][mask]

		for btups in experimental:
			if btups[1] >= args.Tneck:
				phi_diffs = btups[0] - phi_top 
				if len(phi_diffs[phi_diffs>0]) == 0:
					dphi_pos = 1-btups[0]
				else:
					dphi_pos   = np.min(phi_diffs[phi_diffs>0])

				if len(phi_diffs[phi_diffs<0]) == 0:
					dphi_neg = btups[0]
				else:
					dphi_neg   = np.min(np.abs(phi_diffs[phi_diffs<0]))
			else:
				phi_diffs = btups[0] - phi_bottom 
				if len(phi_diffs[phi_diffs>0]) == 0:
					dphi_pos   = 1-btups[0]
				else:
					dphi_pos   = np.min(phi_diffs[phi_diffs>0])

				if len(phi_diffs[phi_diffs<0]) == 0:
					dphi_neg = btups[0]
				else:
					dphi_neg = np.min(np.abs(phi_diffs[phi_diffs<0]))
			weights.append(dphi_pos+dphi_neg)

	elif args.label == "UCST":
		for btups in experimental:
			phi_diffs = btups[0] - experimental[:,0] 
			if len(phi_diffs[phi_diffs>0]) == 0:
				dphi_pos = 1-btups[0]
			else:
				dphi_pos   = np.min(phi_diffs[phi_diffs>0])

			if len(phi_diffs[phi_diffs<0]) == 0:
				dphi_neg = btups[0]
			else:
				dphi_neg   = np.min(np.abs(phi_diffs[phi_diffs<0]))
		weights.append(dphi_pos+dphi_neg)

	weights = np.array(weights)/np.sum(weights)
	return weights

if __name__=="__main__":

	if args.label == "UCST":
		if args.Tucst is None:
			print(f"Label provided is \"UCST\" but no option for --T-UCST is provided. Exiting...", flush=True)
			exit()

	elif args.label == "NECK":
		if args.Tneck is None:
			print(f"Label provided is \"NECK\" but no option for --T-neck is provided. Exiting...", flush=True)
			exit()

	if args.T is None:
		print(f"No temperature range provided to option --T. Exiting...", flush=True)
		exit()

	if args.label not in ["UCST", "NECK", "LOOP"]:
		print(f"Incorrect label provided... Exiting.", flush=True)
		exit()

	start = time.time()

	fig = plt.figure(num=0, figsize=(3,3))
	ax  = plt.axes()
	ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')

	for i in range(len(args.pv)):

		# get the experimental data points
		df           = pd.read_csv(args.csv[i], sep=',', engine="python", names=["phi", "T"])
		experimental = np.array([df["phi"].values, df["T"].values]).T
		
		if args.cw:
			weights = loss_weights(experimental)
		else:
			weights = np.ones(len(experimental))

		# get the spinodal
		def get_spinodal(params):
			
			pv = 1/(1+np.exp(-params[0]))
			pw = 1/(1+np.exp(-params[5]))

			T = np.logspace (np.log10(args.T[0]), np.log10(args.T[1]), int(1e+6) )
			chi = 24/T * (pv * ( fmsa (params[1], params[2], pw, T) * params[1] + \
				(1 - fmsa (params[1], params[2], pw, T) ) * params[2] ) + (1-pv) * params[2] \
				- 1/2 * ( (params[3] + params[4]) ) )

			p1 =  1/(4 * chi * params[-2] * params[-1]) * ( -params[-2] + params[-1] + 2 * chi * params[-2] * params[-1] - np.sqrt(-8 * chi * params[-2] * params[-1] ** 2 + (-params[-2] + params[-1] + 2 * chi * params[-2] * params[-1])**2))
			p2 =  1/(4 * chi * params[-2] * params[-1]) * ( -params[-2] + params[-1] + 2 * chi * params[-2] * params[-1] + np.sqrt(-8 * chi * params[-2] * params[-1] ** 2 + (-params[-2] + params[-1] + 2 * chi * params[-2] * params[-1])**2))

			mask = np.logical_and(np.logical_and(p1>0, p1<1), np.logical_and(p2>0, p2<1))
			p1   = p1[mask]
			p2   = p2[mask]
			T    = T [mask]

			mask = np.logical_and(np.isnan(p1), np.isnan(p2))
			p1   = p1[~mask]
			p2   = p2[~mask]
			T    = T [~mask]

			return p1, p2, T

		# plot the binodal 
		def plot_binodal(params, ax, plot_info):

			pv = 1/(1+np.exp(-params[0]))
			pw = 1/(1+np.exp(-params[5]))

			chi = lambda T: 24/T * (pv * ( fmsa (params[1], params[2], pw, T) * params[1] + \
				(1 - fmsa (params[1], params[2], pw, T) ) * params[2] ) + (1-pv) * params[2] \
				- 1/2 * ( (params[3] + params[4]) ) )

			mu_s = lambda phi, T: np.log(1-phi) + (1 - params[-1]/params[-2]) * phi       + chi(T) * (phi ** 2)     * params[-1]
			mu_p = lambda phi, T: np.log(phi)   + (1 - params[-2]/params[-1]) * (1 - phi) + chi(T) * (1 - phi) ** 2 * params[-2]

			spin_phi_left, spin_phi_right, spin_T = get_spinodal(params)
			print(f"spin_T = {spin_T}")
			if len(spin_T) == 0:
				return 

			else:
				# print(f"Entering binodal computation...")
				arm_left  = []
				arm_right = []
				arm_T     = []

				if args.label == "LOOP":
					for idx, T in enumerate(spin_T[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left [idx*args.skip]/1e+9), np.log10(spin_phi_left[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left[idx*args.skip] and root[1] >= spin_phi_right[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right[idx*args.skip] and root[1] <= spin_phi_left[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				elif args.label == "UCST":
					mask           = spin_T > args.Tucst
					spin_T         = spin_T[mask]
					spin_phi_left  = spin_phi_left[mask]
					spin_phi_right = spin_phi_right[mask]

					print(f"spin_phi_left  = {spin_phi_left}")
					print(f"spin_phi_right = {spin_phi_right}")

					for idx, T in enumerate(spin_T[::args.skip]):
						# print(f"@ T = {T}...", flush=True)

						phi_left  = np.logspace(np.log10(spin_phi_left[idx*args.skip]/1e+9), np.log10(spin_phi_left[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right[idx*args.skip])    , np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left[idx*args.skip] and root[1] >= spin_phi_right[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right[idx*args.skip] and root[1] <= spin_phi_left[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				elif args.label == "NECK":
					mask            = (spin_T <= args.Tneck)
					spin_T_         = np.flip(spin_T[mask])
					spin_phi_left_  = np.flip(spin_phi_left[mask])
					spin_phi_right_ = np.flip(spin_phi_right[mask])

					for idx, T in enumerate(spin_T_[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left_[idx*args.skip]/1e+9), np.log10(spin_phi_left_[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right_[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left_[idx*args.skip] and root[1] >= spin_phi_right_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right_[idx*args.skip] and root[1] <= spin_phi_left_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

					mask            = (spin_T >= args.Tneck)
					spin_T_         = spin_T[mask]
					spin_phi_left_  = spin_phi_left[mask]
					spin_phi_right_ = spin_phi_right[mask]

					for idx, T in enumerate(spin_T_[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left_[idx*args.skip]/1e+9), np.log10(spin_phi_left_[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right_[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left_[idx*args.skip] and root[1] >= spin_phi_right_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found @ T = {T}")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right_[idx*args.skip] and root[1] <= spin_phi_left_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found @ T = {T}")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				arm_left  = np.array(arm_left)
				arm_right = np.array(arm_right)
				arm_T     = np.array(arm_T)

				top_arm = arm_T > args.Tneck
				ax.plot(arm_left[top_arm] , arm_T[top_arm], lw=2, c=plot_info[0], label="_nolabel_")
				ax.plot(arm_right[top_arm], arm_T[top_arm], lw=2, c=plot_info[0], label="_nolabel_")

				bot_arm = arm_T < args.Tneck
				ax.plot(arm_left[bot_arm] , arm_T[bot_arm], lw=2, c=plot_info[0], label="_nolabel_")
				ax.plot(arm_right[bot_arm], arm_T[bot_arm], lw=2, c=plot_info[0], label="_nolabel_")

			return 

		# get the binodal loss function
		def my_binodal_loss(params):
			
			pv = 1/(1+np.exp(-params[0]))
			pw = 1/(1+np.exp(-params[5]))

			max_energy_param = np.max([params[1], params[2], params[3], params[4]])
			params[1] -= max_energy_param
			params[2] -= max_energy_param
			params[3] -= max_energy_param
			params[4] -= max_energy_param

			chi = lambda T: 24/T * (pv * (fmsa(params[1], params[2], pw, T) * params[1] +  \
				(1 - fmsa (params[1], params[2], pw, T)) * params[2]) + (1-pv) * params[2] \
				- 1/2 * (params[3] + params[4]))

			mu_s = lambda phi, T: np.log(1-phi) + (1 - params[-1]/params[-2]) * phi + chi(T)       * (phi ** 2) * params[-1]
			mu_p = lambda phi, T: np.log(phi)   + (1 - params[-2]/params[-1]) * (1 - phi) + chi(T) * (1 - phi) ** 2 * params[-2]

			spin_phi_left, spin_phi_right, spin_T = get_spinodal(params)
			if len(spin_T) == 0:
				loss = 1e+8

			else:
				arm_left  = []
				arm_right = []
				arm_T     = []

				if args.label == "LOOP":
					for idx, T in enumerate(spin_T[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left[idx*args.skip]/1e+9), np.log10(spin_phi_left[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left[idx*args.skip] and root[1] >= spin_phi_right[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right[idx*args.skip] and root[1] <= spin_phi_left[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				elif args.label == "UCST":
					# print(f"Finding binodal loss in UCST!", flush=True)
					mask           = spin_T > args.Tucst
					spin_T         = spin_T[mask]
					spin_phi_left  = spin_phi_left [mask]
					spin_phi_right = spin_phi_right[mask]

					# print(f"spin_phi_left  = {spin_phi_left}")
					# print(f"spin_phi_right = {spin_phi_right}")

					for idx, T in enumerate(spin_T[::args.skip]):

						phi_left  = np.logspace(np.log10(spin_phi_left[idx*args.skip]/1e+9), np.log10(spin_phi_left[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left[idx*args.skip] and root[1] >= spin_phi_right[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right[idx*args.skip] and root[1] <= spin_phi_left[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found for T = {T}!")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				elif args.label == "NECK":
					mask            = (spin_T <= args.Tneck)
					spin_T_         = np.flip(spin_T[mask])
					spin_phi_left_  = np.flip(spin_phi_left[mask])
					spin_phi_right_ = np.flip(spin_phi_right[mask])

					for idx, T in enumerate(spin_T_[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left_[idx*args.skip]/1e+9), np.log10(spin_phi_left_[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right_[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left_[idx*args.skip] and root[1] >= spin_phi_right_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right_[idx*args.skip] and root[1] <= spin_phi_left_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

					mask            = (spin_T >= args.Tneck)
					spin_T_         = spin_T[mask]
					spin_phi_left_  = spin_phi_left[mask]
					spin_phi_right_ = spin_phi_right[mask]

					for idx, T in enumerate(spin_T_[::args.skip]):
						phi_left  = np.logspace(np.log10(spin_phi_left_[idx*args.skip]/1e+9), np.log10(spin_phi_left_[idx*args.skip]), 1000)
						phi_right = np.logspace(np.log10(spin_phi_right_[idx*args.skip]), np.log10(1 - 1e-3), 1000)

						mu_s_left, mu_s_right = mu_s(phi_left, T), mu_s(phi_right, T)
						mu_p_left, mu_p_right = mu_p(phi_left, T), mu_p(phi_right, T)

						mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
						mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)

						s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
						p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]

						def delta_mu(phi):
							eq1 = mu_s(phi[0], T) - mu_s(phi[1], T)
							eq2 = mu_p(phi[0], T) - mu_p(phi[1], T)
							return [eq1, eq2]
						root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)

						if root[0] <= spin_phi_left_[idx*args.skip] and root[1] >= spin_phi_right_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found @ T = {T}")
								arm_left.append(root[0])
								arm_right.append(root[1])
								arm_T.append(T)
							else:
								continue
						
						elif root[0] >= spin_phi_right_[idx*args.skip] and root[1] <= spin_phi_left_[idx*args.skip]:
							if (np.abs(delta_mu(root)) < 1e-6).all():
								# print(f"Found @ T = {T}")
								arm_left.append (root[1])
								arm_right.append(root[0])
								arm_T.append(T)
							else:
								continue
						else:
							continue

				arm_left  = np.array(arm_left)
				arm_right = np.array(arm_right)
				arm_T     = np.array(arm_T)

				if len(arm_T) == 0:
					print(f"There is no temperature in binodal.")
					loss = 1e+8 

				else:
					phis = np.hstack((arm_left, arm_right))
					Ts   = np.hstack((arm_T, arm_T))
					binodal_points = np.array([phis, Ts]).T 
					dsquare = 0
					for idx, point in enumerate(experimental):
						distances = np.linalg.norm(binodal_points - point.reshape(-1), axis=1)
						min_distance = np.min(distances)
						dsquare += min_distance**2
						# delta           = np.linalg.norm(binodal_points[:,0] - point[0]) # find the index at which the phi's match
						# idx_phi_similar = np.argsort(delta) 
						# min_T_distance = np.min(binodal_points[idx_phi_similar[:5],1] - point[1])
						# dsquare += weights[idx] * (min_T_distance**2)

					loss = dsquare 

					if args.label == "LOOP":
						exp_min, exp_max = np.min(experimental[:,1]),   np.max(experimental[:,1])
						bin_min, bin_max = np.min(binodal_points[:,1]), np.max(binodal_points[:,1])
						loss += (exp_min-bin_min)**2 + (exp_max-bin_max)**2
						# print(f"exp_min = {exp_min}, exp_max = {exp_max}")
						# print(f"bin_min = {bin_min}, bin_max = {bin_max}")
				

			return loss

		# get the params 
		params = [args.pv[i], args.emsa[i], args.emsn[i], args.emm[i], args.ess[i], args.pw[i], args.vm[i], args.vs[i]]

		# get the plot info
		plot_info = [args.c[i], args.marker[i]]
		print(f"Plotting binodal...", end=' ', flush=True)
		plot_binodal(params, ax, plot_info)
		ax.scatter(experimental[:,0], experimental[:,1], marker=args.marker[i], s=50, c=args.c[i], edgecolors='k', label="experimental data", zorder=20)
		print(f"done!", flush=True)

	# start plotting everything 
	print(f"Plotting everything...", end=' ', flush=True)
	
	ax.set_ylim(args.ylim[0], args.ylim[1])
	ax.set_xlim(args.xlim[0], args.xlim[1])
	ax.minorticks_on()
	ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
	ax.set_yticks([200, 275, 350, 425, 500])
	# ax.set_xticks(np.array(np.linspace(args.xlim[0], args.xlim[1], 5), dtype=int))
	# ax.set_yticks(np.array(np.linspace(args.ylim[0], args.ylim[1], 5)))
	ax.set_xticklabels([])
	ax.set_yticklabels([])

	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
	print(f"done!", flush=True)

	stop = time.time()
	print(f"Time for computation is {stop - start} seconds.")
