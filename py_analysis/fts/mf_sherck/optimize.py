#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt
import field 
import cma

def linear_B_term(v_pp, v_ss, rho_p, rho_s, N):
	return (v_ss) / (rho_p * N) + v_pp / rho_s

def quadratic_C_term(v_pp, v_ps, v_ss):
	return v_pp * v_ss - v_ps ** 2

def Tmin (B, C):
	return -2 * C/ B

def fmin(A, B, C):
	return A - B**2 / (4 * C)

if __name__=="__main__":
	
	# set up some constants
	rho_p, rho_s, N = 0.1, 1, 50

	# set up the constant term in the quadratic
	A = 1 / (rho_p * rho_s * N)

	# set up the constraint function
	def loss(x):
		v_pp, v_ps, v_ss = x[0], x[1], x[2]
		B = linear_B_term(v_pp, v_ss, rho_p, rho_s, N)
		C = quadratic_C_term(v_pp, v_ps, v_ss)
		Tmin_val = Tmin(B, C)
		fmin_val = fmin(A, B, C)

		# get individual penalties
		Bpenalty    = B**2 if B > 0 else 0
		Cpenalty    = C**2 if C < 0 else 0
		Tminpenalty = Tmin_val**2 if Tmin_val < 0 else 0
		fminpenalty = fmin_val**2 if fmin_val > 0 else 0

		# define the penalty
		penalty = Bpenalty + Cpenalty + Tminpenalty + fminpenalty

		return penalty

	# set up the optimization problem
	x0     = [-0.1, 0.1, 0] # initial guess
	sigma0 = 0.5        # initial step size
	
	# run cma-es
	options = {
		'maxiter': 10000,
		'popsize': 10,
		'verbose': -1,
		'tolfun': 1e-12,       # stop if function value doesn't improve
		'tolfunrel': 1e-9,     # relative function value improvement
		'tolstagnation': 100,  # stop if stuck for 100 generations
	}
	res = cma.fmin(loss, x0, sigma0, options)

	# print the results
	optimized_x0 = res[0] # Best solution
	fval         = res[1] # corresponding loss

	print(f"Optimized x0: {optimized_x0}")
	print(f"Function value: {fval}")
