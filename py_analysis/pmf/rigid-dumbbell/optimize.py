#!/home/satyend/.conda/envs/phase/bin/python

import RigidDumbbell
import numpy as np
import cma

def sigmoid(x):
	return 1.0 / (1.0 + np.exp(-x))

def transform_params(params):
	# params: raw CMA parameters
	# first 3 are epsilons, next 3 are sigmas
	eps_AP = 100.0 * sigmoid(params[0])
	eps_BP = 100.0 * sigmoid(params[1])
	eps_PP = 100.0 * sigmoid(params[2])
	sigma_AP = 100.0 * sigmoid(params[3])
	sigma_BP = 100.0 * sigmoid(params[4])
	sigma_PP = 100.0 * sigmoid(params[5])
	return eps_AP, eps_BP, eps_PP, sigma_AP, sigma_BP, sigma_PP

def objective(raw_params, thermo):
	eps_AP, eps_BP, eps_PP, sigma_AP, sigma_BP, sigma_PP = transform_params(raw_params)
	eps = {"AP": eps_AP, "BP": eps_BP, "PP": eps_PP}
	sig = {"AP": sigma_AP, "BP": sigma_BP, "PP": sigma_PP}

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
	if chi_max <= 4:
		penalty += (4-chi_max)**2
	return penalty

if __name__=="__main__":

	# define an initial guess for epsilons...
	eps_AP = 0.0
	eps_BP = 5.0
	eps_PP = 0.1

	# ... and sigmas
	sigma_AP = 2.0
	sigma_BP = 0.2
	sigma_PP = 0.4

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
	
	# get the temperatures
	# temp_list = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 5, 7.5, 10, 15, 25, 50, 100])
	temp_list = np.array([0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 5, 7.5, 10])

	# get the object down
	thermo = RigidDumbbell.RigidDumbbell(d=0.1, sigmas=sigmas, epsilons=epsilons, temps=temp_list)

	def wrapped_objective(params):
		val = objective(params, thermo)
		return val

	# set up some inputs for cma
	x0 = [-50, -0.569, -4.55404, -0.405, -3.26984, -2.05296]
	sig0 = 1.0
	opts = {
		"popsize": 10,
		"maxiter": 100,
		"ftarget": 1e-12,
		"verb_disp": 1,
		"tolfun": 1e-6,
		"tolx": 1e-6,
		"seed": 42}
	print(f"options = {opts}", flush=True)
	# es = cma.CMAEvolutionStrategy(x0, sig0, opts)
	# es.optimize(wrapped_objective)
	result = cma.fmin(wrapped_objective, x0, sig0, options=opts)
	print(f"Optimized loss = {result[1]}", flush=True)

	# get the best parameters
	best_raw_params = result[0]
	best_params     = transform_params(best_raw_params)

	# FINAL UPDATE
	print(f"Best raw parameters: {best_raw_params}", flush=True)
	print(f"Best parameters: {best_params}", flush=True)
	print(f"Running a final update with these parameters...", flush=True)
	epsilons = {"AP": best_params[0], "BP": best_params[1], "PP": best_params[2]}
	sigmas = {"AP": best_params[3], "BP": best_params[4], "PP": best_params[5]}
	thermo.__reupdate__(sigmas=sigmas, epsilons=epsilons)

	# recompute and plot
	thermo.run()
	thermo.plot()

