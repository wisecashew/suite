import numpy as np
import matplotlib.pyplot as plt
import argparse

# set up the description 
parser = argparse.ArgumentParser(description="Enter name of image file to generate.")
parser.add_argument("-o", metavar="image.png", dest='o', type=str, help="Enter name of file which will be filled with coordinates") 
args   = parser.parse_args()

if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(6,4), nrows=2, ncols=2)
	# We change the fontsize of minor ticks label 
	ax[0][0].tick_params(axis='both', which='major', labelsize=4)
	ax[0][1].tick_params(axis='both', which='major', labelsize=4)
	ax[1][0].tick_params(axis='both', which='major', labelsize=4)
	ax[1][1].tick_params(axis='both', which='major', labelsize=4)

	Nm      = 32
	z       = 26
	emm     = -1
	ems     = -0.541666
	emc     = -0.916666 # -0.583333
	esc     = 0 # -0.416666
	T_list  = [0.1, 1.0, 10.0]
	colors  = ["darkblue", "gold", "darkred"]

	for idx, T in enumerate(T_list):

		chi_ms = 24/T * (ems - 1/2 * emm)
		chi_mc = 24/T * (emc - 1/2 * emm)
		chi_sc = 24/T * (esc) 

		mu_s   = lambda xc: np.log(xc)     + chi_sc * (1 - xc)**2
		mu_c   = lambda xc: np.log(1 - xc) + chi_sc * xc ** 2

		xc_bulk = [0.001, 0.01, 0.1 , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999]
		Ns_solv = []
		Nc_solv = []
		Nm_solv = []
		N_solv  = []

		for xc_b in xc_bulk:
			f = lambda Ntot, xc: T * chi_ms * Nm * (1-xc) + T * chi_mc * Nm * xc + \
			T * Ntot * (1-xc) * np.log(1-xc) + T * Ntot * xc * np.log(xc) - T * Ntot * mu_s(xc_b) * (1-xc) - T * Ntot * mu_c(xc_b) * xc

			NTOT, XC = np.meshgrid(np.linspace(0,z*Nm,1000), np.linspace(0,1,100))
			NTOT = NTOT.flatten()
			XC   = XC.flatten()

			free_energy = f(NTOT, XC) 
			mask        = np.isnan(free_energy)
			free_energy = free_energy[~mask]
			XC          = XC[~mask]
			idx_min     = np.argmin(free_energy)

			Ntot = NTOT[idx_min]
			Ns_solv.append(Ntot*(1-XC[idx_min]))
			Nc_solv.append(Ntot*(XC[idx_min]))
			Nm_solv.append(0)
			N_solv.append(Ntot)

		ax[0][0].plot(xc_bulk, Nm_solv, marker='o', mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", clip_on=False)
		ax[0][1].plot(xc_bulk, Ns_solv, marker='o', mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", clip_on=False)
		ax[1][0].plot(xc_bulk, Nc_solv, marker='o', mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", clip_on=False)
		ax[1][1].plot(xc_bulk, N_solv, marker='o', mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}",  clip_on=False)
		ax[0][0].legend(loc="upper right", prop={'size':5})
		ax[0][1].legend(loc="upper right", prop={'size':5})
		ax[1][0].legend(loc="upper right", prop={'size':5})
		ax[1][1].legend(loc="upper right", prop={'size':5})
		ax[1][0].set_xlabel("Bulk number fraction of cosolvent", fontsize=6)
		ax[1][1].set_xlabel("Bulk number fraction of cosolvent", fontsize=6)
		ax[0][0].set_ylabel("fraction of monomer", fontsize=6) 
		ax[0][1].set_ylabel("fraction of solvent", fontsize=6)
		ax[1][0].set_ylabel("fraction of cosolvent", fontsize=6)
		ax[1][1].set_ylabel("Total particles", fontsize=6)
		ax[0][0].set_xlim(0,1)
		# ax[0][0].set_ylim(0,1)
		ax[0][1].set_xlim(0,1)
		# ax[0][1].set_ylim(0,1)
		ax[1][0].set_xlim(0,1)
		# ax[1][0].set_ylim(0,1)
		ax[1][1].set_xlim(0,1)
		# ax[3].set_ylim(0,1)
	
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
