#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import cma
import sys
sys.path.append("/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Analytical_Solvation/exec/new/build")
import partition_module as pm
import json
from datetime import datetime

class CMAOptimizer:
	def __init__(self, T_arr, target_Nmm, target_Nms_a,
				 Nm=40, Nmm_max=273, alpha=7.333501277066573, 
				 beta=34.26095605400542, dfile="thermo.dump"):
		"""
		Initialize optimizer with target values at a single temperature
		"""
		self.T_arr        = T_arr
		self.target_Nmm   = target_Nmm
		self.target_Nms_a = target_Nms_a
		self.Nm           = Nm
		self.Nmm_max      = Nmm_max
		self.alpha        = alpha
		self.beta         = beta
		self.dfile        = dfile
		
		# Keep track of all evaluations
		self.evaluation_history = []
		return
		
	def objective_function(self, params):
		"""
		Calculate error between model predictions and target values
		params: [EMM_A, EMS_A, EMS_N, pw]
		"""
		total_error = 0
		EMM_A, EMS_A, EMS_N, pw = params
		total_error += ((EMS_N - EMS_A)<0) * np.abs(EMS_N-EMS_A)
		
		ave_Nmm   = []
		ave_Nms_a = []
		
		# Create partition object
		partition = pm.Partition(
			Nm      = self.Nm,
			Nmm_max = self.Nmm_max,
			T       = 1,
			pw      = pw,
			EMM_A   = EMM_A,
			EMM_N   = 0,
			EMS_A   = EMS_A,
			EMS_N   = EMS_N,
			alpha   = self.alpha,
			beta    = self.beta,
			dfile   = self.dfile
		)
		for temp in self.T_arr:
			# set the temperature
			partition.T = temp
			
			# Calculate properties
			partition.get_partition_weights_quick()
			partition.get_partition()
			partition.get_average_Nmm()
			partition.get_average_Nmm_a()
			partition.get_average_Nms_a()
			ave_Nmm.append(partition.ave_Nmm)
			ave_Nms_a.append(partition.ave_Nms_a)
		
		# numpy-ify everything
		ave_Nmm   = np.array(ave_Nmm)
		ave_Nms_a = np.array(ave_Nms_a)

		# Calculate normalized errors
		error_Nmm    = np.linalg.norm((ave_Nmm   - self.target_Nmm)   / self.target_Nmm)   ** 2
		error_Nms_a  = np.linalg.norm((ave_Nms_a - self.target_Nms_a) / self.target_Nms_a) ** 2
		total_error += error_Nmm + error_Nms_a
		
		# Save evaluation to history
		eval_record = {
			"params":    params.tolist(),
			"error":     float(total_error),
			"timestamp": datetime.now().isoformat()
		}
		self.evaluation_history.append(eval_record)
		
		# Save history after each evaluation
		self.save_history()
		
		print(f"Params: {params}")
		print(f"Errors - (error) Nmm: {error_Nmm:.6f}, (error) Nms_a: {error_Nms_a:.6f}")
		print(f"Total Error: {total_error:.6f}\n")
		
		return total_error
	
	def save_history(self, filename='optimization_history.json'):
		"""Save optimization history to file"""
		with open(filename, 'w') as f:
			json.dump(self.evaluation_history, f, indent=2)
	
	def optimize(self, initial_params=None, sigma0=0.5):
		"""
		Run CMA-ES optimization
		
		initial_params: Starting point for optimization
		sigma0: Initial step size
		"""
		if initial_params is None:
			# Default initial parameters
			initial_params = [-1, -1, 0, 0, 0.25]  # Your current values
		
		# Parameter bounds
		bounds = [
			[-100, 0],    # EMM_A
			[-100, 0],    # EMS_A
			[-100, 0],    # EMS_N
			[0.1, 0.99]   # pw
		]
		
		# CMA-ES options
		opts = cma.CMAOptions()
		opts['bounds'] = [[], []]  # [[lb1, lb2, ...], [ub1, ub2, ...]]
		for lb, ub in bounds:
			opts['bounds'][0].append(lb)
			opts['bounds'][1].append(ub)
		
		opts['maxiter'] = 100
		opts['popsize'] = 10   # Smaller population due to expensive evaluation
		opts['verbose'] = -1   # Reduce CMA-ES output
		opts['ftarget'] = 1e-6 # Define the stopping criterion

		# Initialize CMA-ES optimizer
		es = cma.CMAEvolutionStrategy(initial_params, sigma0, opts)
		
		# Run optimization
		print("Starting CMA-ES optimization...\n")
		es.optimize(self.objective_function)
		
		return es.result

if __name__ == "__main__":

	# extract targets
	df = pd.read_csv("contacts.avg.chunk8.dump", sep='\s+', names=["T", "ave_Nmm", "s_Nmm", "ave_Nms", "s_Nms", "ave_Nmsa", "s_Nmsa", "ave_Ns", "s_Ns"], engine="python")

	# Example usage
	optimizer = CMAOptimizer(
		T_arr        = df["T"].values,        # Your target temperature
		target_Nmm   = df["ave_Nmm"].values,  # Your target Nmm
		target_Nms_a = df["ave_Nmsa"].values, # Your target Nmm_a
		Nm           = 40,   # degree of polymerization
		Nmm_max      = 273   # maximum number of contacts
	)
	
	# Start from your current parameters
	initial_guess = [0, -9.38332901, -8.83072076,  0.21210379]
	
	result = optimizer.optimize(initial_params=initial_guess)
	
	print("\nOptimization Results:")
	print(f"Best parameters found: {result.xbest}")
	print(f"Best fitness achieved: {result.fbest}")
	print(f"Number of evaluations: {result.evaluations}")
