#!/home/satyend/.conda/envs/phase/bin/python

import FlexibleDumbbell
import numpy as np

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

	# get the temperatures
	temp_list = np.array([0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2.5, 3, 4, 5, 7.5, 10, 15, 25, 50, 100])

	# get the object down
	print(f"Making the object...", end=' ', flush=True)
	thermo = FlexibleDumbbell.FlexibleDumbbell(k=k_spring, sigmas=sigmas, epsilons=epsilons, temps=temp_list)
	print(f"done.", flush=True)
	thermo._setup_grids()
	thermo._setup_potentials()

	# run the entire calc
	thermo.run()

	# plot the thing
	thermo.plot()

	print(f"Exiting main program!", flush=True)
