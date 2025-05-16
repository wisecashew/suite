#!/home/satyend/.conda/envs/phase/bin/python

import SphericalThermodynamics as st
import numpy as np
import cma

def sigmoid(x):
	return 1.0 / (1.0 + np.exp(-x))

def transform_params(params):
	# params: raw CMA parameters
	# first 3 are epsilons, next 3 are sigmas
	eps_PS   = 10.0 * sigmoid(params[0])
	eps_PP   = 10.0 * sigmoid(params[1])
	eps_SS   = 10.0 * sigmoid(params[2])
	sigma_PS = 2.0 * sigmoid(params[3])
	sigma_PP = 2.0 * sigmoid(params[4])
	sigma_SS = 2.0 * sigmoid(params[5])
	return eps_PS, eps_PP, eps_SS, sigma_PS, sigma_PP, sigma_SS

def objective(raw_params, thermo):
	eps_PS, eps_PP, eps_SS, sigma_PS, sigma_PP, sigma_SS = transform_params(raw_params)
	eps = {"PS": eps_PS,   "PP": eps_PP,   "SS": eps_SS}
	sig = {"PS": sigma_PS, "PP": sigma_PP, "SS": sigma_SS}

	# reupdate everything
	thermo.__reupdate__(sigmas=sig, epsilons=eps)

	# recompute 
	thermo.run()

	# get the low temperature chi
	chi_min = thermo.chi[0]

	# get the maximum value of chi
	chi_max = np.max(thermo.chi)

	# go about defining the penalty
	penalty = 0.0
	if chi_min >= 0:
		penalty += (chi_min)**2
	if chi_max <= 1:
		penalty += (1-chi_max)**2
	return penalty

if __name__=="__main__":

	# define an initial guess for epsilons...
	eps_PS = 0.0
	eps_PP = 0.1
	eps_SS = 0.0

	# ... and sigmas
	sigma_PS = 2.0
	sigma_PP = 0.4
	sigma_SS = 1.0

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
	temp_list = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 5, 7.5, 10, 15, 25, 50, 100])

	# get the object down
	thermo = st.SphericalThermodynamics(d=0.5, sigmas=sigmas, epsilons=epsilons, temps=temp_list)

	def wrapped_objective(params):
		val = objective(params, thermo)
		return val

	# set up some inputs for cma
	x0 = [-50, 0, -50, 50, -1.38629, 0]
	sig0 = 1.0
	opts = {
		"popsize": 10,
		"maxiter": 1000,
		"ftarget": 1e-12,
		"verb_disp": 1,
		"tolfun": 1e-6,
		"tolx": 1e-6,
		"seed": 42}
	print(f"options = {opts}", flush=True)

	# run the fmin
	result = cma.fmin(wrapped_objective, x0, sig0, options=opts)
	print(f"Optimized loss = {result[1]}", flush=True)

	# get the best parameters
	best_raw_params = result[0]
	best_params     = transform_params(best_raw_params)

	# FINAL UPDATE
	print(f"Best raw parameters: {best_raw_params}.", flush=True)
	print(f"Best parameters: {best_params}.", flush=True)
	print(f"Running a final update with these parameters...", flush=True)
	epsilons = {"PS": best_params[0], "PP": best_params[1], "SS": best_params[2]}
	sigmas   = {"PS": best_params[3], "PP": best_params[4], "SS": best_params[5]}
	thermo.__reupdate__(sigmas=sigmas, epsilons=epsilons)

	# recompute and plot
	thermo.run()
	thermo.plot()

