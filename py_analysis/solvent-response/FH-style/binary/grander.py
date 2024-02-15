import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root
import sympy as sym
import argparse 
import warnings
import linecache

EPS=1e-4

parser = argparse.ArgumentParser (description="Plots phase diagrams.")
parser.add_argument ("--vm",   dest='vm',   type=float, action='store', help="Provide vm.", default=-1)
parser.add_argument ("--vs",   dest='vs',   type=float, action='store', help="Provide vs.", default=-1)
parser.add_argument ("--pv",   dest='pv',   type=float, nargs='+', action='store', help="Provide PV.", default=-1)
parser.add_argument ("--emma", dest="emma", type=float, nargs='+', action='store', help="Provide EMMA.", default=-1)
parser.add_argument ("--emmn", dest="emmn", type=float, nargs='+', action='store', help="Provide EMMN.", default=-1)
parser.add_argument ("--emsa", dest="emsa", type=float, nargs='+', action='store', help="Provide EMSA.", default=-1)
parser.add_argument ("--emsn", dest="emsn", type=float, nargs='+', action='store', help="Provide EMSN.", default=0)
parser.add_argument ("--essa", dest="essa", type=float, nargs='+', action='store', help="Provide ESSA.", default=0)
parser.add_argument ("--essn", dest="essn", type=float, nargs='+', action='store', help="Provide ESSN.", default=0)
parser.add_argument ("--pwms", dest="pwms", type=float, nargs='+', action='store', help="Provide PW_MS.", default=0)
parser.add_argument ("--pwmm", dest="pwmm", type=float, nargs='+', action='store', help="Provide PW_MM.", default=0)
parser.add_argument ("--pwss", dest="pwss", type=float, nargs='+', action='store', help="Provide PW_SS.", default=0)
parser.add_argument ("--img",  dest="img",  type=str,   action='store', help="Name of image.", default="bintest.png")
args = parser.parse_args() 

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"Probably a math sqrt error (line {lineno}).\n"

warnings.formatwarning = custom_warning_format

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

zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float64) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float64)
zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float64) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float64)
zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float64) + (1-pw)*np.exp (-1/T * essn, dtype=np.float64)
fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float64)/zmm(emma, emmn, pw, T)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float64)/zms(emsa, emsn, pw, T)
fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float64)/zss(essa, essn, pw, T)

class Phase:

	def __init__ (self, param_list, vm, vs):
		self.EMSA   = param_list[0]
		self.EMSN   = param_list[1]
		self.EMMA   = param_list[2]
		self.EMMN   = param_list[3]
		self.ESSA   = param_list[4]
		self.ESSN   = param_list[5]
		self.PV     = param_list[6]
		self.PWMS   = param_list[7]
		self.PWMM   = param_list[8]
		self.PWSS   = param_list[9]
		self.vm     = vm 
		self.vs     = vs
		return

	def reset_params (self, param_list):
		self.EMSA   = param_list[0]
		self.EMSN   = param_list[1]
		self.EMMA   = param_list[2]
		self.EMMN   = param_list[3]
		self.ESSA   = param_list[4]
		self.ESSN   = param_list[5]
		self.PV     = param_list[6]
		self.PWMS   = param_list[7]
		self.PWMM   = param_list[8]
		self.PWSS   = param_list[9]
		return

	def print_params (self):
		print (f"EMSA = {self.EMSA}, EMSN = {self.EMSN}, EMMA = {self.EMMA}, EMMN = {self.EMMN}, ESSA = {self.ESSA}, ESSN = {self.ESSN}, PV = {self.PV}, PWMS = {self.PWMS}, PWMM = {self.PWMM}, PWSS = {self.PWSS}, vs = {self.vs}, vm = {self.vm}", flush=True)
		return

	def chi (self, T):
		c = 24 * (self.PV * ( (fmsa (self.EMSA, self.EMSN, self.PWMS, T) * self.EMSA + (1 - fmsa (self.EMSA, self.EMSN, self.PWMS, T) ) * self.EMSN) - 1/2 * \
		( (fmma (self.EMMA, self.EMMN, self.PWMM, T) * self.EMMA + (1-fmma (self.EMMA, self.EMMN, self.PWMM, T) ) * self.EMMN) + \
			(fssa (self.ESSA, self.ESSN, self.PWSS, T) * self.ESSA + (1-fssa (self.ESSA, self.ESSN, self.PWSS, T) ) * self.ESSN) ) )
		+ (1-self.PV) * (self.EMSN - 1/2 * (self.EMMN + self.ESSN) ) ) / T

		return c

	def spinodal (self, T):
		p1 =  1/(4 * self.chi(T) * self.vm * self.vs) * ( -self.vm + self.vs + 2 * self.chi(T) * self.vm * self.vs - np.sqrt(-8 * self.chi(T) * self.vm * self.vs ** 2 + (-self.vm + self.vs + 2 * self.chi(T) * self.vm * self.vs)**2))
		p2 =  1/(4 * self.chi(T) * self.vm * self.vs) * ( -self.vm + self.vs + 2 * self.chi(T) * self.vm * self.vs + np.sqrt(-8 * self.chi(T) * self.vm * self.vs ** 2 + (-self.vm + self.vs + 2 * self.chi(T) * self.vm * self.vs)**2))
		
		if isinstance(T, float):
			return [p1, p2, T]
		elif isinstance(T, np.ndarray):
			mask = np.logical_and(np.logical_and(p1>0, p1<1), np.logical_and(p2>0, p2<1))
			p1   = p1[mask]
			p2   = p2[mask]
			T    =  T[mask]
			mask = np.logical_and(np.isnan(p1), np.isnan(p2))
			p1   = p1[~mask]
			p2   = p2[~mask]
			T    = T [~mask]
			return [p1, p2, T]
	
	def spin_condition(self, T, p):
		return -2*self.chi(T) + 1/(p*self.vm) + 1/((1-p)*self.vs)

	def get_critical_info(self):
		phi_  = ( np.sqrt(self.vm*self.vs) - self.vs)/(self.vm - self.vs)
		phi__ = (-np.sqrt(self.vm*self.vs) - self.vs)/(self.vm - self.vs)

		T_up   = np.logspace(1, 3, 5)
		T_down = np.logspace(-2, 2, 5)

		break_for_real = False
		single_T = []

		for Tu in T_up:
			for Td in T_down:
				if phi_ > 0 and phi__>0:
					print(f"phi_ > 0 and phi__>0...")
					self.phi_c = [phi_, phi__]
					self.T_c   = []
					chi_  = (self.vm + 2*(self.vm**(3/2))*(vs**0.5)/(self.vm-self.vs) + 3*self.vs - 2*self.vm*self.vs/(self.vm - self.vs) + 2 * (self.vm**0.5) * self.vs**(3/2)/(self.vm-self.vs) + 2 * self.vs**2/(self.vm-self.vs))/(2*self.vm*self.vs)
					chi__ = (self.vm - 2*(self.vm**(3/2))*(vs**0.5)/(self.vm-self.vs) + 3*self.vs - 2*self.vm*self.vs/(self.vm - self.vs) - 2 * (self.vm**0.5) * self.vs**(3/2)/(self.vm-self.vs) + 2 * self.vs**2/(self.vm-self.vs))/(2*self.vm*self.vs)
					self.chi_c = np.array([chi_, chi__])
					root1      = fsolve(self.spin_condition, Tu,  args=(phi_))
					root2      = fsolve(self.spin_condition, Td, args=(phi__))
					if (np.abs(self.spin_condition(root1[0], phi_)).all() < 1e-6):
						self.T_c.append(root1[0])
					if (np.abs(self.spin_condition(root2[0], phi__)).all() < 1e-6):
						self.T_c.append(root2[0])
					self.T_c   = np.array(self.T_c)
							
				elif phi_ > 0 and phi__ < 0:
					print(f"phi_ > 0 and phi__<0...")
					self.phi_c = np.array([phi_, phi_])
					self.T_c   = []
					chi_       = (self.vm + 2*(self.vm**(3/2))*(vs**0.5)/(self.vm-self.vs) + 3*self.vs - 2*self.vm*self.vs/(self.vm - self.vs) + 2 * (self.vm**0.5) * self.vs**(3/2)/(self.vm-self.vs) + 2 * self.vs**2/(self.vm-self.vs))/(2*self.vm*self.vs)
					self.chi_c = np.array([chi_, chi_])
					root1      = fsolve(self.spin_condition, Tu, args=(phi_))
					root2      = fsolve(self.spin_condition, Td, args=(phi_))
					if (np.abs(self.spin_condition(root1[0], phi_)).all() < 1e-6):
						self.T_c.append(root1[0])
					if (np.abs(self.spin_condition(root2[0], phi_)).all() < 1e-6):
						self.T_c.append(root2[0])
					self.T_c   = np.array(self.T_c)

				elif phi_ < 0 and phi__ > 0:
					print(f"phi_ < 0 and phi__ > 0...")
					self.phi_c = np.array([phi__, phi__])
					self.T_c   = []
					chi__      = (self.vm - 2*(self.vm**(3/2))*(vs**0.5)/(self.vm-self.vs) + 3*self.vs - 2*self.vm*self.vs/(self.vm - self.vs) - 2 * (self.vm**0.5) * self.vs**(3/2)/(self.vm-self.vs) + 2 * self.vs**2/(self.vm-self.vs))/(2*self.vm*self.vs)
					self.chi_c = np.array([chi__, chi__])
					root1      = fsolve(self.spin_condition, Tu , args=(phi__))
					root2      = fsolve(self.spin_condition, Td,  args=(phi__))
					if (np.abs(self.spin_condition(root1[0], phi__)).all() < 1e-6):
						self.T_c.append(root1[0])
					if (np.abs(self.spin_condition(root2[0], phi__)).all() < 1e-6):
						self.T_c.append(root2[0])
					self.T_c   = np.array(self.T_c)

				self.T_c, keep = delete_close_elements(self.T_c)
				if len(self.T_c) == 2:
					break_for_real = True 
					print("yo!")
					break
				
				elif len(self.T_c) == 1:
					single_T.append(self.T_c[0])

				else:
					pass

			if break_for_real:
				break

		if not break_for_real and len(single_T) != 0:
			self.T_c = np.array([single_T[-1]])

		print(self.T_c)
		print(self.chi_c)
		self.T_c, keep = delete_close_elements(self.T_c)
		if len(keep) != 0:
			self.phi_c = self.phi_c[keep]
			self.chi_c = self.chi_c[keep]
			

		return

	def setup(self):
		phi_p1, phi_p2, T_ = sym.symbols('phi_p1 phi_p2 T_')
		zmm   = self.PWMM*sym.exp (-1/T_ * self.EMMA) + (1-self.PWMM)*sym.exp (-1/T_ * self.EMMN)
		zms   = self.PWMS*sym.exp (-1/T_ * self.EMSA) + (1-self.PWMS)*sym.exp (-1/T_ * self.EMSN)
		zss   = self.PWSS*sym.exp (-1/T_ * self.ESSA) + (1-self.PWSS)*sym.exp (-1/T_ * self.ESSN)
		fmma  = self.PWMM*sym.exp (-1/T_ * self.EMMA)/zmm
		fmsa  = self.PWMS*sym.exp (-1/T_ * self.EMSA)/zms
		fssa  = self.PWSS*sym.exp (-1/T_ * self.ESSA)/zss

		chi  = 24/T_ * (self.PV * ( (fmsa * self.EMSA + (1 - fmsa)  * self.EMSN) - 1/2 * \
		( (fmma * self.EMMA + (1-fmma) * self.EMMN) + \
		(fssa * self.ESSA + (1-fssa) * self.ESSN) ) ) \
		+ (1-self.PV) * (self.EMSN - 1/2 * (self.EMMN + self.ESSN) ) )

		mu_s1 = sym.log(1-phi_p1) + (1 - self.vs/self.vm) * phi_p1 + chi * (phi_p1 ** 2) * self.vs
		mu_s2 = sym.log(1-phi_p2) + (1 - self.vs/self.vm) * phi_p2 + chi * (phi_p2 ** 2) * self.vs

		mu_p1 = sym.log(phi_p1) + (1 - self.vm/self.vs) * (1 - phi_p1) + chi * (1 - phi_p1) ** 2 * self.vm
		mu_p2 = sym.log(phi_p2) + (1 - self.vm/self.vs) * (1 - phi_p2) + chi * (1 - phi_p2) ** 2 * self.vm

		# lambdify the mu calculations
		L_mu_s = sym.lambdify([phi_p1, T_], mu_s1)
		L_mu_p = sym.lambdify([phi_p1, T_], mu_p1)

		self.mu_s = lambda phi_p, temp: L_mu_s(phi_p, temp)
		self.mu_p = lambda phi_p, temp: L_mu_p(phi_p, temp)

		delta_mu_s = mu_s1 - mu_s2 
		delta_mu_p = mu_p1 - mu_p2 

		# lambdify the functions
		L_delta_mu_s = sym.lambdify([phi_p1, phi_p2, T_], delta_mu_s)
		L_delta_mu_p = sym.lambdify([phi_p1, phi_p2, T_], delta_mu_p)
		
		self.delta_mu_s = lambda p1, p2, temp: L_delta_mu_s(p1, p2, temp)
		self.delta_mu_p = lambda p1, p2, temp: L_delta_mu_p(p1, p2, temp)

		# do the differentiation
		d_delta_mu_s_dpp1 = sym.diff(delta_mu_s, phi_p1)
		d_delta_mu_s_dpp2 = sym.diff(delta_mu_s, phi_p2)
		d_delta_mu_s_T    = sym.diff(delta_mu_s, T_)

		# lambdify the functions
		L_d_delta_mu_s_dpp1 = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_s_dpp1)
		L_d_delta_mu_s_dpp2 = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_s_dpp2)
		L_d_delta_mu_s_T    = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_s_T)

		# get the functions
		self.d_delta_mu_s_dpp1 = lambda p1, p2, temp: L_d_delta_mu_s_dpp1(p1, p2, temp)
		self.d_delta_mu_s_dpp2 = lambda p1, p2, temp: L_d_delta_mu_s_dpp2(p1, p2, temp)
		self.d_delta_mu_s_T    = lambda p1, p2, temp: L_d_delta_mu_s_T(p1, p2, temp)

		# do the differentiation 
		d_delta_mu_p_dpp1 = sym.diff(delta_mu_p, phi_p1)
		d_delta_mu_p_dpp2 = sym.diff(delta_mu_p, phi_p2)
		d_delta_mu_p_T    = sym.diff(delta_mu_p, T_) 

		# lambdify the functions
		L_d_delta_mu_p_dpp1 = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_p_dpp1)
		L_d_delta_mu_p_dpp2 = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_p_dpp2)
		L_d_delta_mu_p_T    = sym.lambdify([phi_p1, phi_p2, T_], d_delta_mu_p_T)

		# get the functions 
		self.d_delta_mu_p_dpp1 = lambda p1, p2, temp: L_d_delta_mu_p_dpp1(p1, p2, temp)
		self.d_delta_mu_p_dpp2 = lambda p1, p2, temp: L_d_delta_mu_p_dpp2(p1, p2, temp)
		self.d_delta_mu_p_T    = lambda p1, p2, temp: L_d_delta_mu_p_T(p1, p2, temp)

	def get_nbrhd_solution(self, T_c, phi_c):
		delta_phi_p2 = np.logspace(-2,-6, 100)
		for delta_pp2 in delta_phi_p2:
			def diff_mu(p):
				eq1 = self.delta_mu_s(p[0], phi_c+delta_pp2, p[1])
				eq2 = self.delta_mu_p(p[0], phi_c+delta_pp2, p[1])
				return [eq1, eq2]

			root1 = fsolve(diff_mu, [phi_c-delta_pp2, T_c+0.01], xtol=1e-6)
			root2 = fsolve(diff_mu, [phi_c-delta_pp2, T_c-0.01], xtol=1e-6)
			if (np.abs(diff_mu(root1)) < 1e-6).all():
				return root1[0], phi_c + delta_pp2, root1[1]

			elif (np.abs(diff_mu(root2)) < 1e-6).all():
				return root2[0], phi_c + delta_pp2, root2[1]
			
			else:
				continue 

		return None

	def calc_perturbations(self, bin_arm_1, bin_arm_2, T_list, delta_T):
		a1 = self.d_delta_mu_s_dpp1(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])
		b1 = self.d_delta_mu_s_dpp2(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])
		c1 = -delta_T * self.d_delta_mu_s_T(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])

		a2 = self.d_delta_mu_p_dpp1(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])
		b2 = self.d_delta_mu_p_dpp2(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])
		c2 = -delta_T * self.d_delta_mu_p_T(bin_arm_1[-1], bin_arm_2[-1], T_list[-1])

		D  = np.linalg.det(np.array([[a1, b1], [a2, b2]]))
		Dx = np.linalg.det(np.array([[c1, b1], [c2, b2]]))
		Dy = np.linalg.det(np.array([[a1, c1], [a2, c2]]))

		return Dx/D, Dy/D

	def get_binodal(self, T_arm, bin_arm_left, bin_arm_right):

		T_idx     = index_max_distance(T_arm) 
		T_list    = [T_arm[T_idx-1], T_arm[T_idx]]
		bin_arm_1 = [bin_arm_left [T_idx-1], bin_arm_left[T_idx]]
		bin_arm_2 = [bin_arm_right[T_idx-1], bin_arm_right[T_idx]]
		delta_T   = T_list[-1] - T_list[-2] 
		iterx = 0
		scale = False

		while iterx < 100000 and T_list[-1] < T_arm[T_idx+1] and delta_T > 1e-12:

			if scale:
				delta_T /= 1.1
			else:
				delta_T = T_list[-1] - T_list[-2] 
			if iterx % 1000 == 0:
				print(f"iterx: {iterx}, delta_T = {delta_T}")
			delta_phi_p1, delta_phi_p2 = self.calc_perturbations(bin_arm_1, bin_arm_2, T_list, delta_T)

			def delta_mu(phi):
				eq1 = self.delta_mu_s(phi[0], phi[1], T_list[-1]+delta_T)
				eq2 = self.delta_mu_p(phi[0], phi[1], T_list[-1]+delta_T)
				return [eq1, eq2]
			
			root = fsolve(delta_mu, [bin_arm_1[-1]+delta_phi_p1, bin_arm_2[-1]+delta_phi_p2])
			if (np.array(delta_mu(root)) < 1e-6).all():
				scale = False
				bin_arm_1.append(root[0])
				bin_arm_2.append(root[1])
				T_list.append(T_list[-1]+delta_T)
			else:
				scale = True
			iterx += 1
		

		return bin_arm_1, bin_arm_2, T_list

	def grid_search(self):
		if len(self.T_c) == 1:
			T_range = np.logspace(np.log10(self.T_c[0]/1000), np.log10(self.T_c[0]),10000)
		elif len(self.T_c) == 2:
			phi_range = self.spinodal(T_range[0] - 0.0001)
			if ~np.isnan(phi_range[0]) and ~np.isnan(phi_range[1]): # isinstance(phi_range[0], float):
				T_range = np.linspace(phase.T_c[0]/1000, phase.T_c[1], 10000)
			elif isinstance(phi_range[0], np.ndarray):
				T_range = np.logspace(np.log10(self.T_c[0]/1000), np.log10(self.T_c[1]),10000)
		
		# get p boundaries at each temperature
		arm_left  = []
		arm_right = []
		arm_T     = []
		for Tidx, T in enumerate(T_range[1:-1]):
			if Tidx % 1000 == 0:
				print(f"idx = {Tidx}", flush=True)
			phi_range  = self.spinodal(T)
			phi_left   = np.logspace(np.log10(phi_range[0]/1000000), np.log10(phi_range[0]), 1000)
			phi_right  = np.logspace(np.log10(phi_range[1]), 0, 1000)
			mu_s_left, mu_s_right = self.mu_s(phi_left,  T), self.mu_s(phi_right, T)
			mu_p_left, mu_p_right = self.mu_p(phi_left,  T), self.mu_p(phi_right, T)
			mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
			mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)
			s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
			p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]
			def delta_mu(phi):
				eq1 = self.delta_mu_s(phi[0], phi[1], T)
				eq2 = self.delta_mu_p(phi[0], phi[1], T)
				return [eq1, eq2]
			root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)
			if root[0] < phi_range[0] and root[1] > phi_range[1]:
				if (np.abs(delta_mu(root)) < 1e-6).all():
					arm_left.append (root[0])
					arm_right.append(root[1])
					arm_T.append(T)
				else:
					continue 
			elif root[0] > phi_range[1] and root[1] < phi_range[0]:
				if (np.abs(delta_mu(root)) < 1e-6).all():
					arm_left.append (root[1])
					arm_right.append(root[0])
					arm_T.append(T)
				else:
					continue
			else:
				continue

		return arm_left, arm_right, arm_T

	def grid_search_neo(self, T_arm, bin_arm_left, bin_arm_right):
		# get p boundaries at each temperature
		arm_left  = []
		arm_right = []
		arm_T     = []
		for idx, T in enumerate(T_arm):
			if idx % 1000 == 0:
				print(f"idx = {idx}", flush=True)
			phi_left  = np.logspace(np.log10(bin_arm_left[idx]/1e+9), np.log10(bin_arm_left[idx]), 1000)
			phi_right = np.logspace(np.log10(bin_arm_right[idx]), np.log10(1 - (1-(bin_arm_right[idx]))/1e+9), 1000)
			mu_s_left, mu_s_right = self.mu_s(phi_left, T), self.mu_s(phi_right, T)
			mu_p_left, mu_p_right = self.mu_p(phi_left, T), self.mu_p(phi_right, T)
			mu_s_left_idx, mu_s_right_idx, min_mus_dist = find_closest_indices_vectorized(mu_s_left, mu_s_right)
			mu_p_left_idx, mu_p_right_idx, min_mup_dist = find_closest_indices_vectorized(mu_p_left, mu_p_right)
			s_phi_left, s_phi_right = phi_left[mu_s_left_idx], phi_right[mu_s_right_idx]
			p_phi_left, p_phi_right = phi_left[mu_p_left_idx], phi_right[mu_p_right_idx]
			def delta_mu(phi):
				eq1 = self.delta_mu_s(phi[0], phi[1], T)
				eq2 = self.delta_mu_p(phi[0], phi[1], T)
				return [eq1, eq2]
			root = fsolve(delta_mu, [(s_phi_left+p_phi_left)/2, (s_phi_right+p_phi_right)/2], xtol=1e-6)
			if root[0] < bin_arm_left[idx] and root[1] > bin_arm_right[idx]:
				if (np.abs(delta_mu(root)) < 1e-6).all():
					arm_left.append(root[0])
					arm_right.append(root[1])
					arm_T.append(T)
				else:
					# print(f"No sol for T = {T}... min: dmus_dist = {min_mus_dist}, dmup_dist = {min_mup_dist} ")
					# print(f"Spin left = {bin_arm_left[idx]}, spin right = {bin_arm_right[idx]}")
					# print(f"Left guess = {(s_phi_left+p_phi_left)/2}, right guess = {(s_phi_right+p_phi_right)/2}")
					continue
			elif root[0] > bin_arm_right[idx] and root[1] < bin_arm_left[idx]:
				if (np.abs(delta_mu(root)) < 1e-6).all():
					arm_left.append (root[1])
					arm_right.append(root[0])
					arm_T.append(T)
				else:
					# print(f"No sol for T = {T}... min: dmus_dist = {min_mus_dist}, dmup_dist = {min_mup_dist} ")
					# print(f"Spin left = {bin_arm_left[idx]}, spin right = {bin_arm_right[idx]}")
					# print(f"Left guess = {(s_phi_left+p_phi_left)/2}, right guess = {(s_phi_right+p_phi_right)/2}")
					continue
			else:
				# print(f"Wack roots...")
				# print(f"No sol for T = {T}...")
				# print(f"Spin left = {bin_arm_left[idx]}, spin right = {bin_arm_right[idx]}")
				# print(f"Left guess = {(s_phi_left+p_phi_left)/2}, right guess = {(s_phi_right+p_phi_right)/2}")
				continue

		return arm_left, arm_right, arm_T

if __name__=="__main__":

	fig = plt.figure(num=0, figsize=(3,3))
	ax  = plt.axes()
	vm  = args.vm
	vs  = args.vs

	params = [args.emsa[0], args.emsn[0], args.emma[0], args.emmn[0], args.essa[0], args.essn[0], args.pv[0], args.pwms[0], args.pwmm[0], args.pwss[0]]

	phase = Phase(params, vm, vs) 
	phase.print_params()
	phase.setup()
	phase.get_critical_info()

	print(f"Critical temperatures = {phase.T_c}")
	print(f"Critical phi          = {phase.phi_c}")

	phase.T_c = np.sort(phase.T_c)

	if len(phase.T_c) == 1:
		T_range = np.logspace(np.log10(phase.T_c[0]/1000), np.log10(phase.T_c[0]),10000)
	elif len(phase.T_c) == 2:
		phi_range = phase.spinodal(phase.T_c[0] - EPS)
		if ~np.isnan(phi_range[0]) and ~np.isnan(phi_range[1]):
			T_range = np.linspace(phase.T_c[0]/1000, phase.T_c[1], 10000)
		else:
			T_range = np.logspace(np.log10(phase.T_c[0]/1000), np.log10(phase.T_c[1]),10000)

	print(f"Final T_range = {T_range}")
	arms    = phase.spinodal(T_range)

	print(f"arms[0] = {arms[0]}", flush=True)
	print(f"arms[1] = {arms[1]}", flush=True)
	
	ax.scatter(phase.phi_c, phase.T_c, c='darkred', s=1, zorder=10)
	ax.scatter(arms[0], arms[2], s=0.5, c='coral', zorder=1)
	ax.scatter(arms[1], arms[2], s=0.5, c='steelblue', zorder=1)

	
	print(f"Start looking for binodals...")
	arms[0] = np.flip(arms[0])
	mask = arms[0] > arms[0][0]/10
	arms[0] = arms[0][mask]
	arms[1] = np.flip(arms[1])
	arms[1] = arms[1][mask]
	arms[2] = np.flip(arms[2])
	arms[2] = arms[2][mask]

	'''
	arm_left, arm_right, T_list  = phase.grid_search_neo(arms[2], arms[0], arms[1])
	ax.scatter(arm_left,  T_list, s=0.5, c='slategray', zorder=1)
	ax.scatter(arm_right, T_list, s=0.5, c='black',     zorder=1)

	bin_arm_left, bin_arm_right, T_list = phase.get_binodal(T_list, arm_left, arm_right)
	ax.scatter(bin_arm_left,  T_list, s=0.5, c='lavender', zorder=1)
	ax.scatter(bin_arm_right, T_list, s=0.5, c='pink',     zorder=1)
	'''
	# df = pd.read_csv("PEG-water.csv", sep=',', engine="python", names=["phi", "T"])
	# ax.scatter(df["phi"].values, df["T"].values, marker='^', s=5, c='gold', edgecolors='k')

	# ax.set_ylim(420, 540)
	# ax.set_yticks(np.arange(420, 540, 20))
	# ax.set_xlim(-0.05, 0.5)

	fig.savefig(args.img, dpi=1200, bbox_inches="tight")
	
