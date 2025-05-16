#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numba
from numba import njit, prange

# define the different factories when determining potentials
# ===========================================================
# this is for the square well
def square_well_factory(eps, sigma):
	def potential(r):
		return np.where(r < sigma, eps, 0.0)
	return potential

# this is for the gaussian
def gaussian_factory(eps, sigma):
	if eps <= 0 or sigma <= 0:
		print(f"You cannot have negative epsilon or sigma when writing a Gaussian potential.", flush=True)
		print(f"Exiting...", flush=True)
		exit()
	def potential(r):
		gauss = eps/(4*np.pi*sigma**2)**(3/2) * np.exp(-(r**2)/(4*sigma**2))
		return gauss
	return potential
# ===========================================================

@njit(inline='always')
def _v_gauss(r, eps, sig):
	pref = eps / ((4.0*np.pi*sig**2)**1.5)
	return pref * np.exp(-(r*r)/(4.0*sig*sig))

@njit(parallel=True, fastmath=True)
def _pmf_PS_flexible_numba(R_vals, theta_grid, phi_grid, 
				sin_theta,
				dtheta, dphi,
				axis_x, axis_y, axis_z,
				sigma_AP, eps_AP,
				sigma_BP, eps_BP,
				beta, k_spring, l_nodes, l_weights):
	"""
	Vectorised + parallel PMF.  Returns W_PS(R) for every R_vals.
	"""
	nR   = R_vals.size
	n_th = theta_grid.shape[0]
	n_ph = theta_grid.shape[1]
	n_l  = l_nodes.size

	# harmonic prefactor (beta*k/2pi)^(3/2)
	harm_pref = ((beta * k_spring) / (2.0*np.pi))**(3/2)

	# this is the pmf
	W = np.empty(nR, dtype=np.float64)

	for i in prange(nR):                                # parallel over radial shells
		R  = R_vals[i]
		xP = R                                          # P lies on +x axis
		Z  = 0.0                                        # rotational partition function

		for j in range(n_l):                            # quadrature over the angles
			# at this particular bond value... 
			l     = l_nodes[j]
			w_lag = l_weights[j]

			spring_w = harm_pref * np.exp(-0.5 * beta * k_spring * l * l)
			w_tot    = w_lag * spring_w * l * l

			for ith in range(n_th):                         # pure Python loops compiled by numba
				for iph in range(n_ph):

					# A & B positions
					ax = axis_x[ith, iph]
					ay = axis_y[ith, iph]
					az = axis_z[ith, iph]

					# rapid because everything is scalar in these loops
					xA = -0.5 * l * ax;  yA = -0.5 * l * ay;  zA = -0.5 * l * az
					xB =  0.5 * l * ax;  yB =  0.5 * l * ay;  zB =  0.5 * l * az

					# distances to P
					rA = np.sqrt((xP-xA)**2 + yA**2 + zA**2)
					rB = np.sqrt((xP-xB)**2 + yB**2 + zB**2)

					U  = _v_gauss(rA, eps_AP, sigma_AP) + \
						_v_gauss(rB, eps_BP, sigma_BP)

					Z += w_tot * np.exp(-beta*U) * sin_theta[ith, iph]

		Z   *= dtheta*dphi / (4.0*np.pi)
		W[i] = -(1.0/beta)*np.log(Z)

	# make sure there are no arbitrary constant hanging around
	W = W - np.min(W)
	return W

class FlexibleDumbbell:
	def __init__(self,
				k: float,
				sigmas: dict,
				epsilons: dict,
				temps: list,
				R_max: float = 5.0,
				N_R:  int = 500,
				N_th: int = 90,
				N_ph: int = 180,
				N_l:  int = 50):
		# set up the properties of the object
		self.k          = k
		self.sigmas     = sigmas
		self.epsilons   = epsilons
		self.temps      = temps

		# set up the discretization information
		self.R_max = R_max
		self.N_R   = N_R
		self.N_th  = N_th
		self.N_ph  = N_ph
		self.N_l   = N_l
		return

	def __reupdate__(self, sigmas, epsilons, k_spring):
		self.sigmas   = sigmas
		self.epsilons = epsilons
		self.k        = k_spring
		self._setup_potentials()
		return

	def _setup_potentials(self):
		self.potentials = {
			"AP": gaussian_factory(self.epsilons["AP"], self.sigmas["AP"]),
			"BP": gaussian_factory(self.epsilons["BP"], self.sigmas["BP"]),
			"PP": gaussian_factory(self.epsilons["PP"], self.sigmas["PP"])
		}
		return

	def _setup_length_grid(self, beta: float):
		# set up the rigid, old model
		if self.k is None:
			self.l_vals = np.array([self.d]) 
			self.w_l    = np.array([1.0])
		
		else:
			self.l_vals = np.linspace(0.0, self.l_cut, self.N_l, endpoint=False)
			rho = (self.l_vals**2)*np.exp(-beta*self.k*self.l_vals**2/2)
			self.w_l = rho/np.sum(rho) # normalize
		
		return

	def _setup_grids(self):
		sigma_min   = min(self.sigmas.values())
		self.R_min  = 0.001 * sigma_min
		self.R_vals = np.logspace(np.log10(self.R_min), np.log10(self.R_max), self.N_R)
		self.theta  = np.linspace(0, np.pi, self.N_th)
		self.phi    = np.linspace(0, 2*np.pi, self.N_ph, endpoint=False)
		self.dtheta = self.theta[1] - self.theta[0]
		self.dphi   = self.phi[1] - self.phi[0]
		self.theta_grid, self.phi_grid = np.meshgrid(self.theta, self.phi, indexing='ij')
		self.sin_theta = np.sin(self.theta_grid)

		# get the axes
		self.axis_x = np.sin(self.theta_grid)*np.cos(self.phi_grid)
		self.axis_y = np.sin(self.theta_grid)*np.sin(self.phi_grid)
		self.axis_z = np.cos(self.theta_grid)

		# get the l grid
		self.x_nodes, self.x_weights = np.polynomial.laguerre.laggauss(self.N_l)
		return

	def compute_pmf_PS(self, T: float) -> np.ndarray:
		beta = 1.0 / T
		scale = np.sqrt(2/(beta * self.k))
		l_nodes   = scale * np.sqrt(self.x_nodes)
		l_weights = 0.5 * scale / np.sqrt(self.x_nodes) * self.x_weights
		W_PS_3D = _pmf_PS_flexible_numba(self.R_vals, self.theta_grid, self.phi_grid,
										self.sin_theta, self.dtheta, self.dphi,
										self.axis_x, self.axis_y, self.axis_z,
										self.sigmas["AP"], self.epsilons["AP"],
										self.sigmas["BP"], self.epsilons["BP"],
										beta, self.k, l_nodes, l_weights)
		return W_PS_3D

	def compute_pmf_PP(self, T:float) -> np.ndarray:
		W_PP_3D = self.potentials["PP"](self.R_vals)
		return W_PP_3D

	def compute_B2_PS(self, W_PS_3D: np.ndarray, T: float) -> float:
		beta = 1.0 / T
		f    = np.exp(-beta * W_PS_3D) - 1
		B2   = -2*np.pi * np.trapz(f * self.R_vals**2, self.R_vals)
		return B2

	def compute_B2_PP(self, W_PP_3D, T) -> np.ndarray:
		beta = 1.0 / T
		f    = np.exp(-beta * W_PP_3D) - 1
		B2   = -2*np.pi * np.trapz(f * self.R_vals**2, self.R_vals)
		return B2

	def run(self):
		# set up the grids
		self._setup_grids()
		self._setup_potentials()

		# initialize containers
		b2_PS    = []
		b2_PP    = []
		W_PS_all = []
		W_PP_all = []

		# print(f"\tRunning the temperature loop for potentials of mean force and B2...", end=' ', flush=True)
		for idx,T in enumerate(self.temps):
			# get the potentials of mean force
			W_PS_3D = self.compute_pmf_PS(T)
			W_PP_3D = self.compute_pmf_PP(T)

			# append all comptutation to their respective lists 
			b2_PS.append(self.compute_B2_PS(W_PS_3D, T))
			b2_PP.append(self.compute_B2_PP(W_PP_3D, T))
			W_PS_all.append(W_PS_3D)
			W_PP_all.append(W_PP_3D)
		# print(f"done.", flush=True)

		# numpy-ify everything
		b2_PS = np.array(b2_PS)
		b2_PP = np.array(b2_PP)

		# now get chi_PS
		chi = b2_PS - 0.5 * b2_PP
		self.W_PS_all = W_PS_all
		self.W_PP_all = W_PP_all
		self.b2_PS    = b2_PS
		self.b2_PP    = b2_PP
		self.chi      = chi
		return

	def plot(self):
		# set up the figures
		print(f"\tPlotting the thermodynamic summary...", end=' ', flush=True)
		fig, axs = plt.subplots(2, 2, figsize=(8, 8))
		for idx, T in enumerate(self.temps):
			if idx % 5 == 0:
				axs[0,0].plot(self.R_vals[::5], self.W_PS_all[idx][::5], marker='o', lw=0.5, mec='k', markersize=6, label=f"T = {T}")
				axs[0,1].plot(self.R_vals[::5], self.W_PP_all[idx][::5], marker='o', lw=0.5, mec='k', markersize=6, label=f"T = {T}")

		# get the plots down
		# for W_PS
		axs[0,0].legend(loc="best")
		axs[0,0].set_xlabel("Radial distance $r$")
		axs[0,0].set_ylabel("$W_{PS}(r)$")
		axs[0,0].grid(visible=True)

		# for W_PP
		axs[0,1].legend(loc="best")
		axs[0,1].set_xlabel("Radial distance $r$")
		axs[0,1].set_ylabel("$W_{PP}(r)$")
		axs[0,1].grid(visible=True)

		# for b2_PS
		axs[1,0].plot(self.temps, self.b2_PS, c="darkred", marker="^", markersize=6, lw=1, ls='--', mec='k')
		axs[1,0].set_xscale("log")
		axs[1,0].set_xlabel("Temperature")
		axs[1,0].set_ylabel("$B_{2,PS}(T)$")
		axs[1,0].grid(visible=True)

		# for chi
		axs[1,1].plot(self.temps, self.chi, c="salmon", marker="o", markersize=6, lw=1, ls='--', mec='k')
		axs[1,1].set_xscale("log")
		axs[1,1].set_xlabel("Temperature")
		axs[1,1].set_ylabel("$\\chi(T)$")
		axs[1,1].grid(visible=True)
		fig.tight_layout()
		fig.savefig("thermodynamics_summary.png", dpi=1200)
		print(f"done.", flush=True)
		return

