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
parser.add_argument ("--T",      dest='T',     type=float, nargs='+', action='store', help="Provide a temperature range to plot thing in.", default=None)
parser.add_argument ("--T-loop", dest='Tloop', type=float, action='store', help="Provide a temperature to bisect a loop.", default=None)
parser.add_argument ("--T-UCST", dest='Tucst', type=float, action='store', help="Provide a temperature above which to hunt for a binodal.", default=None)
parser.add_argument ("--T-neck", dest='Tneck', type=float, action='store', help="Provide a temperature above and below which we will hunt for a binodal.", default=None)
parser.add_argument ("--custom-weights", dest='cw',   action='store_true', help="Enter option to draw custom weights.", default=False)
parser.add_argument ("--optimize",       dest='opt',  action='store_true', help="Enter option to run the optimization.", default=False)
parser.add_argument ("--maxiter",        dest='mi',   type=int, action='store', help="Provide maximal number of iterations.", default=250)
parser.add_argument ("--skip",           dest='skip', type=int, action='store', help="Provide length of jump while looping through points for the binodal.", default=1000)
parser.add_argument ("--inp",    dest='inp',   type=str,   action='store', help='Name of input parameter file.')
parser.add_argument ("--oup",    dest='oup',   type=str,   action='store', help='Name of output parameter file.', default=None)
parser.add_argument ("--csv",    dest='csv',   type=str,   action='store', help='Name of csv file.')
parser.add_argument ("--img",    dest="img",   type=str,   action='store', help="Name of image.", default="bintest.png")
parser.add_argument ("--label",  dest='label', type=str,   action='store', help="Prove the label for the diagram.")
parser.add_argument ("--Tbot",   dest='Tbot',  type=float, action='store', help="Provide a temperature below which you will not search for a binodal.", default=None)
parser.add_argument ("--Ttop",   dest='Ttop',  type=float, action='store', help="Provide a temperature above which you will not search for a binodal.", default=None)
parser.add_argument ("--L1",  dest="L1",  type=float, action='store', help="Provide L1 regularization.", default=0.001)
parser.add_argument ("--L2",  dest="L2",  type=float, action='store', help="Provide L2 regularization.", default=0.001)
parser.add_argument ("--sig", dest="sig", type=float, action='store', help="Provide sigma.", default=0.01)
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

def read_input(input_file):
	f = open (input_file, 'r')
	for line in f:
		L = line.strip().split()
		if L[0] == "pv":
			pv = float(L[2])
		elif L[0] == "pw":
			pw = float(L[2])
		elif L[0] == "emsa":
			emsa = float(L[2])
		elif L[0] == "emsn":
			emsn = float(L[2])
		elif L[0] == "emma":
			emma = float(L[2])
		elif L[0] == "emmn":
			emmn = float(L[2])
		elif L[0] == "essa":
			essa = float(L[2])
		elif L[0] == "essn":
			essn = float(L[2])
		elif L[0] == "vm":
			vm = float(L[2])
		elif L[0] == "vs":
			vs = float(L[2])
		else:
			print(f"Bad input! Exiting...")
			f.close()
			exit()
	f.close()

	params = [pv, pw, emsa, emsn, emma, emmn, essa, essn, vm, vs]
	return params

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

zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma) + (1-pw)*np.exp (-1/T * emmn)
zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa) + (1-pw)*np.exp (-1/T * emsn)
zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa) + (1-pw)*np.exp (-1/T * essn)
fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma)/zmm(emma, emmn, pw, T)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa)/zms(emsa, emsn, pw, T)
fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa)/zss(essa, essn, pw, T)

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

	# get the L1 and L2 hyperparameters
	L1 = args.L1 
	L2 = args.L2

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

	# get the experimental data points
	df           = pd.read_csv(args.csv, sep=',', engine="python", names=["phi", "T"])
	experimental = np.array([df["phi"].values, df["T"].values]).T
	if args.cw:
		weights = loss_weights(experimental)
	else:
		weights = np.ones(len(experimental))

	print(f"weights = {weights}")

	# get the spinodal
	def get_spinodal(params):
		
		# tweak the params
		pv   = 1/(1+np.exp(-params[0]))
		pw   = 1/(1+np.exp(-params[1]))
		emsa = params[2]
		emsn = params[3]
		emma = params[4]
		emmn = params[5]
		essa = params[6]
		essn = params[7]
		vm   = params[8]
		vs   = params[9]

		# temperature range
		T = np.logspace (np.log10(args.T[0]), np.log10(args.T[1]), int(1e+6) )

		# get the FHP chi
		chi = 24/T * (pv * ( (fmsa (emsa, emsn, pw, T) * emsa + (1 - fmsa (emsa, emsn, pw, T) ) * emsn) - 1/2 * \
		( (fmma (emma, emmn, pw, T) * emma + (1-fmma (emma, emmn, pw, T) ) * emmn) + \
		(fssa (essa, essn, pw, T) * essa + (1-fssa (essa, essn, pw, T) ) * essn) ) ) \
		+ (1-pv) * (emsn - 1/2 * (emmn + essn) ) )

		p1 =  1/(4 * chi * vm * vs) * ( -vm + vs + 2 * chi * vm * vs - np.sqrt(-8 * chi * vm * vs ** 2 + (-vm + vs + 2 * chi * vm * vs)**2))
		p2 =  1/(4 * chi * vm * vs) * ( -vm + vs + 2 * chi * vm * vs + np.sqrt(-8 * chi * vm * vs ** 2 + (-vm + vs + 2 * chi * vm * vs)**2))

		mask = np.logical_and(np.logical_and(p1>0, p1<1), np.logical_and(p2>0, p2<1))
		p1   = p1[mask]
		p2   = p2[mask]
		T    = T[mask]

		mask = np.logical_and(np.isnan(p1), np.isnan(p2))
		p1   = p1[~mask]
		p2   = p2[~mask]
		T    = T [~mask]

		return p1, p2, T

	# plot the spinodal 
	def plot_spinodal(params, ax, plot_info):

		p1, p2, T = get_spinodal(params)
		ax.plot(p1, T, lw=0.5, c=plot_info[1], label=plot_info[0], marker=plot_info[3], markersize=0.5)
		ax.plot(p2, T, lw=0.5, c=plot_info[2], label="_nolabel_",  marker=plot_info[3], markersize=0.5)

		return 

	# plot the binodal 
	def plot_binodal(params, ax, plot_info):

		pv = 1/(1+np.exp(-params[0]))
		pw = 1/(1+np.exp(-params[1]))
		emsa = params[2]
		emsn = params[3]
		emma = params[4]
		emmn = params[5]
		essa = params[6]
		essn = params[7]
		vm   = params[8]
		vs   = params[9]

		chi = lambda T: 24/T * (pv * ( (fmsa (emsa, emsn, pw, T) * emsa + (1 - fmsa (emsa, emsn, pw, T) ) * emsn) - 1/2 * \
		( (fmma (emma, emmn, pw, T) * emma + (1-fmma (emma, emmn, pw, T) ) * emmn) + \
		(fssa (essa, essn, pw, T) * essa + (1-fssa (essa, essn, pw, T) ) * essn) ) ) \
		+ (1-pv) * (emsn - 1/2 * (emmn + essn) ) )

		mu_s = lambda phi, T: np.log(1-phi) + (1 - vs/vm) * phi + chi(T) * (phi ** 2) * vs
		mu_p = lambda phi, T: np.log(phi)   + (1 - vm/vs) * (1 - phi) + chi(T) * (1 - phi) ** 2 * vm

		spin_phi_left, spin_phi_right, spin_T = get_spinodal(params)
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
				mask = spin_T > args.Tucst
				spin_T         = spin_T[mask]
				spin_phi_left  = spin_phi_left[mask]
				spin_phi_right = spin_phi_right[mask]

				for idx, T in enumerate(spin_T[::args.skip]):
					# print(f"@ T = {T}...", flush=True)

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
			# print(f"arm_left = {arm_left}")
			# print(f"arm_right = {arm_right}")
			# print(f"arm_T = {arm_T}")
			ax.plot(arm_left , arm_T, lw=0.5, c=plot_info[1], label=plot_info[0], marker=plot_info[3], markersize=0.5)
			ax.plot(arm_right, arm_T, lw=0.5, c=plot_info[2], label="_nolabel_",  marker=plot_info[3], markersize=0.5)

		return 

	# get the spinodal loss function
	def my_spinodal_loss(params):

		pv = 1/(1+np.exp(-params[0]))
		pw = 1/(1+np.exp(-params[1]))
		emsa = params[2]
		emsn = params[3]
		emma = params[4]
		emmn = params[5]
		essa = params[6]
		essn = params[7]
		vm   = params[8]
		vs   = params[9]

		# print(f"pv = {pv}, pw = {pw}")

		T = np.logspace (np.log10(args.T[0]), np.log10(args.T[1]), int(1e+6) )
		chi = 24/T * (pv * ( (fmsa (emsa, emsn, pw, T) * emsa + (1 - fmsa (emsa, emsn, pw, T) ) * emsn) - 1/2 * \
		( (fmma (emma, emmn, pw, T) * emma + (1-fmma (emma, emmn, pw, T) ) * emmn) + \
		(fssa (essa, essn, pw, T) * essa + (1-fssa (essa, essn, pw, T) ) * essn) ) ) \
		+ (1-pv) * (emsn - 1/2 * (emmn + essn) ) )

		p1 =  1/(4 * chi * vm * vs) * ( -vm + vs + 2 * chi * vm * vs - np.sqrt(-8 * chi * vm * vs ** 2 + (-vm + vs + 2 * chi * vm * vs)**2))
		p2 =  1/(4 * chi * vm * vs) * ( -vm + vs + 2 * chi * vm * vs + np.sqrt(-8 * chi * vm * vs ** 2 + (-vm + vs + 2 * chi * vm * vs)**2))

		mask = np.logical_and(np.logical_and(p1>0, p1<1), np.logical_and(p2>0, p2<1))
		p1   = p1[mask]
		p2   = p2[mask]
		T    = T[mask]

		mask = np.logical_and(np.isnan(p1), np.isnan(p2))
		p1   = p1[~mask]
		p2   = p2[~mask]
		T    = T [~mask]

		if len(T) == 0:
			loss = 1e+16

		else:
			phis = np.hstack((p1, p2))
			Ts   = np.hstack((T , T))
			spinodal_points = np.array([phis, Ts]).T

			dsquare = 0
			for idx, point in enumerate(experimental):
				distances = np.linalg.norm(spinodal_points - point.reshape(-1), axis=1)
				min_distance = np.min(distances)
				dsquare += min_distance**2

			loss = dsquare 

		return loss

	# get the binodal loss function
	def my_binodal_loss(params):
		
		pv = 1/(1+np.exp(-params[0]))
		pw = 1/(1+np.exp(-params[1]))
		emsa = params[2]
		emsn = params[3]
		emma = params[4]
		emmn = params[5]
		essa = params[6]
		essn = params[7]
		vm   = params[8]
		vs   = params[9]

		energy = np.array([emsa, emsn, emma, emmn, essa, essn])

		chi = lambda T: 24/T * (pv * ( (fmsa (emsa, emsn, pw, T) * emsa + (1 - fmsa (emsa, emsn, pw, T) ) * emsn) - 1/2 * \
		( (fmma (emma, emmn, pw, T) * emma + (1-fmma (emma, emmn, pw, T) ) * emmn) + \
		(fssa (essa, essn, pw, T) * essa + (1-fssa (essa, essn, pw, T) ) * essn) ) ) \
		+ (1-pv) * (emsn - 1/2 * (emmn + essn) ) )

		mu_s = lambda phi, T: np.log(1-phi) + (1 - vs/vm) * phi       + chi(T) * (phi ** 2)     * vs
		mu_p = lambda phi, T: np.log(phi)   + (1 - vm/vs) * (1 - phi) + chi(T) * (1 - phi) ** 2 * vm

		spin_phi_left, spin_phi_right, spin_T = get_spinodal(params)

		if len(spin_T) == 0:
			loss = 1e+16

		else:
			arm_left  = []
			arm_right = []
			arm_T     = []
			loss      = 0

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
				
				# get an appropriate loop error
				arm_left  = np.array(arm_left)
				arm_right = np.array(arm_right)
				arm_T     = np.array(arm_T)

				top_half  = arm_T >  args.Tloop
				bot_half  = arm_T <= args.Tloop
				arm_T_top     = arm_T    [top_half]
				arm_left_top  = arm_left [top_half]
				arm_right_top = arm_right[top_half]
				arm_T_bot     = arm_T    [bot_half]
				arm_left_bot  = arm_left [bot_half]
				arm_right_bot = arm_right[bot_half]

				exp_top = experimental[experimental[:,1] > args.Tloop]
				exp_bot = experimental[experimental[:,1] <= args.Tloop]
				_loss   = 0

				if len(arm_T_top) > 0:
					for pt in exp_top:
						left_dist  = np.min(np.abs(arm_left_top  - pt[0]))
						right_dist = np.min(np.abs(arm_right_top - pt[0]))
						if left_dist > right_dist:
							min_phi_idx = np.argmin(np.abs(arm_right_top - pt[0]))
							_loss += (arm_T_top[min_phi_idx] - pt[1]) ** 2
						else:
							min_phi_idx = np.argmin(np.abs(arm_left_top - pt[0]))
							_loss += (arm_T_top[min_phi_idx] - pt[1]) ** 2
				
				if len(arm_T_bot) > 0:
					for pt in exp_bot:
						left_dist  = np.min(np.abs(arm_left_bot  - pt[0]))
						right_dist = np.min(np.abs(arm_right_bot - pt[0]))
						if left_dist > right_dist:
							min_phi_idx = np.argmin(np.abs(arm_right_bot - pt[0]))
							_loss += (arm_T_bot[min_phi_idx] - pt[1]) ** 2
						else:
							min_phi_idx = np.argmin(np.abs(arm_left_bot - pt[0]))
							_loss += (arm_T_bot[min_phi_idx] - pt[1]) ** 2

			elif args.label == "UCST":
				mask = spin_T > args.Tucst
				spin_T         = spin_T[mask]
				spin_phi_left  = spin_phi_left[mask]
				spin_phi_right = spin_phi_right[mask]

				for idx, T in enumerate(spin_T[::args.skip]):
					# print(f"@ T = {T}...", flush=True)

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

				# get an appropriate loop error
				arm_left  = np.array(arm_left)
				arm_right = np.array(arm_right)
				arm_T     = np.array(arm_T)

				_loss = 0
				if len(arm_T) > 0:
					for pt in experimental:
						left_dist  = np.min(np.abs(arm_left  - pt[0]))
						right_dist = np.min(np.abs(arm_right - pt[0]))
						if left_dist > right_dist:
							min_phi_idx = np.argmin(np.abs(arm_right - pt[0]))
							_loss += (arm_T[min_phi_idx] - pt[1]) ** 2
						else:
							min_phi_idx = np.argmin(np.abs(arm_left - pt[0]))
							_loss += (arm_T[min_phi_idx] - pt[1]) ** 2

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

				arm_left_bot  = np.array(arm_left)
				arm_right_bot = np.array(arm_right)
				arm_T_bot     = np.array(arm_T)
				exp_bot       = experimental[experimental[:,1] <= args.Tneck]

				# print(f"arm_left_bot  = {arm_left_bot}")
				# print(f"arm_right_bot = {arm_right_bot}")
				# print(f"arm_T_bot     = {arm_T_bot}")
				# print(f"exp_bot       = {exp_bot}")

				_loss = 0
				if len(arm_T_bot) > 0:
					for pt in exp_bot:
						left_dist  = np.min(np.abs(arm_left_bot  - pt[0]))
						right_dist = np.min(np.abs(arm_right_bot - pt[0]))
						if left_dist > right_dist:
							min_phi_idx = np.argmin(np.abs(arm_right_bot - pt[0]))
							# print(f"for pt = {pt}, closest point = ({arm_right_bot[min_phi_idx]},{arm_T_bot[min_phi_idx]}), deltaT = {arm_T_bot[min_phi_idx] - pt[1]}")
							_loss += (arm_T_bot[min_phi_idx] - pt[1]) ** 2
						else:
							min_phi_idx = np.argmin(np.abs(arm_left_bot - pt[0]))
							# print(f"for pt = {pt}, closest point = ({arm_left_bot[min_phi_idx]},{arm_T_bot[min_phi_idx]}), deltaT = {arm_T_bot[min_phi_idx] - pt[1]}")
							_loss += (arm_T_bot[min_phi_idx] - pt[1]) ** 2

				mask            = (spin_T > args.Tneck)
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

				top_half = arm_T > args.Tneck
				arm_left_top  = arm_left[top_half]
				arm_right_top = arm_right[top_half]
				arm_T_top     = arm_T[top_half]
				exp_top       = experimental[experimental[:,1] > args.Tneck]

				if len(arm_T_top) > 0:
					for pt in exp_top:
						left_dist  = np.min(np.abs(arm_left_top  - pt[0]))
						right_dist = np.min(np.abs(arm_right_top - pt[0]))
						if left_dist > right_dist:
							min_phi_idx = np.argmin(np.abs(arm_right_top - pt[0]))
							_loss += (arm_T_top[min_phi_idx] - pt[1]) ** 2
						else:
							min_phi_idx = np.argmin(np.abs(arm_left_top - pt[0]))
							_loss += (arm_T_top[min_phi_idx] - pt[1]) ** 2


			if len(arm_T) == 0:
				loss = 1e+16

			else:
				pass 
				'''
				phis = np.hstack((arm_left, arm_right))
				Ts   = np.hstack((arm_T, arm_T))
				binodal_points = np.array([phis, Ts]).T 

				dsquare = 0
				for idx, point in enumerate(experimental):
					distances = np.linalg.norm(binodal_points - point.reshape(-1), axis=1)
					min_distance = np.min(distances)
					dsquare += weights[idx] * (min_distance**2)
				loss = dsquare 
				'''

				loss    += _loss 
				reg_loss = L1 * np.sum(np.abs(energy)) + L2 * np.sum(energy ** 2)
				loss    +=  reg_loss 
				print(f"reg_loss = {reg_loss}.")
				if emsa > emsn:
					loss += emsa-emsn
				if emma > emmn:
					loss += emma-emmn
				if essa > essn:
					loss += essa-essn

		return loss

	# get the params 
	params = read_input(args.inp) # [args.pv, args.pw, args.emsa, args.emsn, args.emma, args.emmn, args.essa, args.essn, args.vm, args.vs]

	# plot the initial spinodal
	plot_info = ["spinodal (initial)", "crimson", "crimson", 'o']
	print(f"Plotting spinodal...", end=' ', flush=True)
	plot_spinodal(params, ax, plot_info)
	print(f"done!", flush=True)

	plot_info = ["binodal (initial)", "indianred", "indianred", 'o']
	print(f"Plotting binodal...", end=' ', flush=True)
	plot_binodal(params, ax, plot_info)
	print(f"done!", flush=True)

	L_spin = my_spinodal_loss(params)
	print(f"initial spinodal loss is {L_spin}.", flush=True)

	L_bin  = my_binodal_loss(params)
	print(f"initial binodal loss is {L_bin}.", flush=True)

	if args.opt:
		# define the initial guess and std dev
		initial_guess = params 
		print(f"Dimensionality of search space is {len(initial_guess)}.", flush=True)
		initial_sigma = args.sig

		# define the options for the CMA-ES algorithm
		options = {
			'popsize': 10,   # population size of candidate solutions
			'tolfun' : 1e-6, # termination criterion: tolerance in functions  
			'tolx'   : 1e-6, # tolerance criterion: tolerance in x-values during the iterations
			'maxiter': args.mi, # maximum number of iterations to perform
			'verb_disp': 100, # display info every 'verb_disp' iterations
		}

		# run the CMA-ES algorithm
		print(f"Running CMA-ES...", flush=True)
		best_solution = cma.fmin(my_binodal_loss, initial_guess, initial_sigma, options=options)
		print(f"done!", flush=True)

		if args.oup is None:
			pass
		else:
			f = open(args.oup, 'w')
			for i in range(len(best_solution[0])):
				if i == 0:
					f.write(f"pv = {best_solution[0][i]}\n")
				elif i == 1:
					f.write(f"pw = {best_solution[0][i]}\n")
				elif i == 2:
					f.write(f"emsa = {best_solution[0][i]}\n")
				elif i == 3:
					f.write(f"emsn = {best_solution[0][i]}\n")
				elif i == 4:
					f.write(f"emma = {best_solution[0][i]}\n")
				elif i == 5:
					f.write(f"emmn = {best_solution[0][i]}\n")
				elif i == 6:
					f.write(f"essa = {best_solution[0][i]}\n")
				elif i == 7:
					f.write(f"essn = {best_solution[0][i]}\n")
				elif i == 8:
					f.write(f"vm = {best_solution[0][i]}\n")
				elif i == 9:
					f.write(f"vs = {best_solution[0][i]}\n")
			f.close()

		print(f"The best solution: {best_solution}", flush=True)
		plot_info = ["spinodal (final)", "coral", "coral", '^']
		plot_spinodal(best_solution[0], ax, plot_info)
		plot_info = ["binodal (final)", "gold", "gold", '^']
		plot_binodal(best_solution[0], ax, plot_info)	

	# start plotting everything 
	print(f"Plotting everything...", end=' ', flush=True)
	ax.scatter(experimental[:,0], experimental[:,1], marker='^', s=8, c='gold', edgecolors='k', label="experimental data")
	ax.set_ylim(args.T[0], args.T[1])
	ax.set_xlim(0, 1)
	ax.minorticks_on()
	ax.grid(axis='both')
	ax.legend(prop={"size": 6}, loc='upper right')
	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
	print(f"done!", flush=True)




	stop = time.time()
	print(f"Time for computation is {stop - start} seconds.")
