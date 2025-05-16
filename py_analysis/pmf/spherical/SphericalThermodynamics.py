#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

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

class SphericalThermodynamics:
	def __init__(self,
				d: float,
				sigmas: dict,
				epsilons: dict,
				temps: list,
				R_max: float = 5.0,
				N_R:  int = 500,
				N_th: int = 90,
				N_ph: int = 180):
		# set up the properties of the object
		self.d          = d
		self.sigmas     = sigmas
		self.epsilons   = epsilons
		self.temps      = temps

		# set up the discretization information
		self.R_max = R_max
		self.N_R   = N_R
		self.N_th  = N_th
		self.N_ph  = N_ph
		return

	def __reupdate__(self, sigmas, epsilons):
		self.sigmas   = sigmas
		self.epsilons = epsilons
		self._setup_potentials()
		return

	def _setup_potentials(self):
		self.potentials = {
			"PS": square_well_factory(self.epsilons["PS"], self.sigmas["PS"]),
			"PP": square_well_factory(self.epsilons["PP"], self.sigmas["PP"]),
			"SS": square_well_factory(self.epsilons["SS"], self.sigmas["SS"])
		}
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
		return

	def compute_pmf_PS(self, T: float) -> np.ndarray:
		W_PS_3D = self.potentials["PS"](self.R_vals)
		return W_PS_3D

	def compute_pmf_PP(self, T:float) -> np.ndarray:
		W_PP_3D = self.potentials["PP"](self.R_vals)
		return W_PP_3D

	def compute_pmf_SS(self, T:float) -> np.ndarray:
		W_SS_3D = self.potentials["SS"](self.R_vals)
		return W_SS_3D

	def compute_B2_PS(self, W_PS_3D: np.ndarray, T: float) -> float:
		beta = 1.0 / T
		f    = np.exp(-beta * W_PS_3D) - 1
		B2   = -2*np.pi * np.trapz(f * self.R_vals**2, self.R_vals)
		return B2

	def compute_B2_PP(self, W_PP_3D, T) -> float:
		beta = 1.0 / T
		f    = np.exp(-beta * W_PP_3D) - 1
		B2   = -2*np.pi * np.trapz(f * self.R_vals**2, self.R_vals)
		return B2

	def compute_B2_SS(self, W_SS_3D, T) -> float:
		beta = 1.0 / T
		f    = np.exp(-beta * W_SS_3D) - 1
		B2   = -2*np.pi * np.trapz(f * self.R_vals**2, self.R_vals)
		return B2

	def run(self):
		# set up the grids
		self._setup_grids()
		self._setup_potentials()

		# initialize containers
		b2_PS    = []
		b2_PP    = []
		b2_SS    = []
		W_PS_all = []
		W_PP_all = []
		W_SS_all = []

		for idx,T in enumerate(self.temps):
			# get the potentials of mean force
			W_PS_3D = self.compute_pmf_PS(T)
			W_PP_3D = self.compute_pmf_PP(T)
			W_SS_3D = self.compute_pmf_SS(T)

			# append all comptutation to their respective lists 
			b2_PS.append(self.compute_B2_PS(W_PS_3D, T))
			b2_PP.append(self.compute_B2_PP(W_PP_3D, T))
			b2_SS.append(self.compute_B2_SS(W_SS_3D, T))
			W_PS_all.append(W_PS_3D)
			W_PP_all.append(W_PP_3D)
			W_SS_all.append(W_SS_3D)

		# numpy-ify everything
		b2_PS = np.array(b2_PS)
		b2_PP = np.array(b2_PP)
		b2_SS = np.array(b2_SS)

		# now get chi_PS
		chi = b2_PS - 0.5 * (b2_PP + b2_SS)
		self.W_PS_all = W_PS_all
		self.W_PP_all = W_PP_all
		self.b2_PS    = b2_PS
		self.b2_PP    = b2_PP
		self.b2_SS    = b2_SS
		self.chi      = chi
		return

	def plot(self):
		# set up the figures
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
		return

if __name__ == "__main__":

	# define an initial guess for epsilons...
	# PS, PP, SS
	eps_PS = 1e-12
	eps_PP = 0.20838575957175878
	eps_SS = 0

	# ... and sigmas
	# PS, PP, SS
	sigma_PS = 2.0
	sigma_PP = 0.5687661866419248
	sigma_SS = 1

	# define the dictionaries
	epsilons = {
		"PS": eps_PS,
		"PP": eps_PP,
		"SS": eps_SS
	}
	sigmas = {
		"PS": sigma_PS,
		"PP": sigma_PP,
		"SS": sigma_SS
	}
	# the dictionaries have been defined
	
	# get the temperatures
	temp_list = np.array([0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 4, 5, 7.5, 10, 15, 25, 50, 100])

	# get the object down
	thermo = SphericalThermodynamics(d=0.5, sigmas=sigmas, epsilons=epsilons, temps=temp_list)

	# define the temperature list
	thermo.run()

	# plot the final quantities
	thermo.plot()
