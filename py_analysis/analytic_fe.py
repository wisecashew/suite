#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt

def f (n_, t_, k_, b_, p_, x_):
	t1 = 2 * np.log (x_) 
	t2 = - 6 * (np.sin(p_*np.pi/(2*n_))**2) / (np.pi * (b_ ** 2)) * (x_ ** 2)
	t3 = np.log (24 * np.sqrt(6/np.pi) * (np.sin(np.pi*p_/(2*n_))/b_) ** 3 )

	return -k_*t_*(t1+t2+t3)



if __name__=="__main__":

	N = 21
	b = 1
	T = 2/3
	k = 1

	
	fig = plt.figure (figsize=(5, 5))
	ax  = plt.axes   ()

	p   = 1
	x   = np.linspace (0.001, 20, 1000)

	F   = f (N, T, k, b, p, x)

	ax.plot (x, F, markersize=0, lw=1, color="coral", label="exact")
	ax.legend ()
	ax.set_xlim (0, 20)
	ax.set_ylim (0, 6)
	ax.set_xticks (np.linspace (0, 20, 5))
	ax.set_yticks (np.linspace (0, 6, 4))

	fig.savefig ("exact_fe", dpi=1200, bbox_inches="tight")
