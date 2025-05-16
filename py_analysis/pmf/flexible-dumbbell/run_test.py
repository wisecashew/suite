#!/home/satyend/.conda/envs/phase/bin/python

import time
import numpy as np
import numba
from numba import njit, prange
import matplotlib 
matplotlib.use("agg")
import matplotlib.pyplot as plt
import FlexibleDumbbell

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

	# pre-compute the orientation unit-vectors once
	# axis_x = np.sin(theta_grid)*np.cos(phi_grid)
	# axis_y = np.sin(theta_grid)*np.sin(phi_grid)
	# axis_z = np.cos(theta_grid)

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


if __name__=="__main__":

	# define an initial guess for epsilons...
	eps_AP = 1e-12
	eps_BP = 7.2255725435313
	eps_PP = 0.20838575957175878

	# ... and sigmas
	sigma_AP = 2.0
	sigma_BP = 0.18310171612772858
	sigma_PP = 0.5687661866419248

	# define the dictionaries
	epsilons = {
		"AP": eps_AP,
		"BP": eps_BP,
		"PP": eps_PP
	}
	sigmas = {
		'AP': sigma_AP,
		'BP': sigma_BP,
		'PP': sigma_PP
	}
	# the dictionaries have been defined
	k_spring = 0.01
	n_l      = 50
	x_nodes, x_weights = np.polynomial.laguerre.laggauss(n_l)

	# get the temperatures
	temp_list = np.array([0.0001, 0.001, 0.01, 0.05, 0.1]) # , 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 4, 5, 7.5, 10, 15, 25, 50, 100])

	# get the object down
	thermo = FlexibleDumbbell.FlexibleDumbbell(k=k_spring, sigmas=sigmas, epsilons=epsilons, temps=temp_list)
	thermo._setup_grids()
	thermo._setup_potentials()

	# get the axes
	axis_x = np.sin(thermo.theta_grid)*np.cos(thermo.phi_grid)
	axis_y = np.sin(thermo.theta_grid)*np.sin(thermo.phi_grid)
	axis_z = np.cos(thermo.theta_grid)

	# plot the PMF
	fig, axs = plt.subplots(1, 1, figsize=(6, 6))
	
	# set the timer
	print(f"Starting time...", flush=True)
	start_loop = time.time()

	# run the PMF
	for idx,temp in enumerate(temp_list):
		start     = time.time()
		beta      = 1.0 / thermo.temps[idx]
		scale     = np.sqrt(2/(beta * k_spring))
		l_nodes   = scale * np.sqrt(x_nodes)					# these are the x-values
		l_weights = 0.5 * scale / np.sqrt(x_nodes) * x_weights	# these are the associated weights

		W = _pmf_PS_flexible_numba(thermo.R_vals, thermo.theta_grid, thermo.phi_grid,
						thermo.sin_theta,
						thermo.dtheta, thermo.dphi,
						axis_x, axis_y, axis_z,
						thermo.sigmas["AP"], thermo.epsilons["AP"],
						thermo.sigmas["BP"], thermo.epsilons["BP"],
						beta, k_spring, l_nodes, l_weights)
		if idx % 1 == 0:
			axs.plot(thermo.R_vals, W, marker="^", markersize=6, lw=1, ls='--', mec='k', mew=0.1, label=f"$T={temp}$")
		stop = time.time()
		print(f"\tTime for numba computation is {stop-start:.2f} seconds.", flush=True)

	stop_loop = time.time()
	axs.legend(loc='best', fontsize=8)
	axs.set_xlabel("Radial distance $r$")
	axs.set_ylabel("$W_{PS}(r)$ (numba)")
	axs.grid(visible=True)
	print(f"Time for the entire loop is {stop_loop-start_loop:.2f} seconds.", flush=True)
	fig.tight_layout()
	fig.savefig("comparison.png", dpi=1200)
