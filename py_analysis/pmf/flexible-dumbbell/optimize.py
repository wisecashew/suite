#!/home/satyend/.conda/envs/phase/bin/python

import FlexibleDumbbell
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("agg")
import cma

def sigmoid(x):
	return 1.0 / (1.0 + np.exp(-x))

def unsigmoidify(raw_x0):
	raw_x0[0] = np.log(raw_x0[0]/(20-raw_x0[0]))
	raw_x0[1] = np.log(raw_x0[1]/(20-raw_x0[1]))
	raw_x0[2] = np.log(raw_x0[2]/(20-raw_x0[2]))
	raw_x0[3] = np.log(raw_x0[3]/(5-raw_x0[3]))
	raw_x0[4] = np.log(raw_x0[4]/(5-raw_x0[4]))
	raw_x0[5] = np.log(raw_x0[5]/(5-raw_x0[5]))
	raw_x0[6] = np.log(raw_x0[6]/(10-raw_x0[6]))
	return raw_x0

def transform_params(params):
	# params: raw CMA parameters
	# first 3 are epsilons, next 3 are sigmas, final is the harmonic spring
	eps_AP = 20.0 * sigmoid(params[0])
	eps_BP = 20.0 * sigmoid(params[1])
	eps_PP = 20.0 * sigmoid(params[2])
	sigma_AP = 5.0 * sigmoid(params[3])
	sigma_BP = 5.0 * sigmoid(params[4])
	sigma_PP = 5.0 * sigmoid(params[5])
	k_spring = 10.0 * sigmoid(params[6])
	return eps_AP, eps_BP, eps_PP, sigma_AP, sigma_BP, sigma_PP, k_spring

def objective(raw_params, thermo):
	eps_AP, eps_BP, eps_PP, sigma_AP, sigma_BP, sigma_PP, k_spring = transform_params(raw_params)
	eps = {"AP": eps_AP, "BP": eps_BP, "PP": eps_PP}
	sig = {"AP": sigma_AP, "BP": sigma_BP, "PP": sigma_PP}

	# reupdate everything
	thermo.__reupdate__(sigmas=sig, epsilons=eps, k_spring=k_spring)

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
	
	print(f"chi_min, chi_max = {chi_min:.2f}, {chi_max:.2f}; penalty = {penalty}", flush=True)

	return penalty

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

	# ... and a spring constant
	k_spring = 0.01

	# get the temperatures
	temp_list = np.array([0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 5, 7.5, 10, 15, 25])

	# get the object down
	thermo = FlexibleDumbbell.FlexibleDumbbell(k=k_spring, sigmas=sigmas, epsilons=epsilons, temps=temp_list)
	thermo._setup_grids()
	thermo._setup_potentials()

	def wrapped_objective(params):
		val = objective(params, thermo)
		return val

	# initial params
	raw_x0 = [eps_AP, eps_BP, eps_PP, sigma_AP, sigma_BP, sigma_PP, k_spring]
	x0     = unsigmoidify(raw_x0)

	# set up some inputs for cma
	sig0 = 1.0
	opts = {
		"popsize": 10,
		"maxiter": 10,
		"ftarget": 1e-12,
		"verb_disp": 1,
		"tolfun": 1e-6,
		"tolx": 1e-6,
		"seed": 42}
	print(f"options = {opts}", flush=True)
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
	sigmas   = {"AP": best_params[3], "BP": best_params[4], "PP": best_params[5]}
	k_spring = best_params[6]
	thermo.__reupdate__(sigmas=sigmas, epsilons=epsilons, k_spring=k_spring)

	# recompute and plot
	thermo.run()
	thermo.plot()

