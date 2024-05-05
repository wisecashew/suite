import numpy as np
import matplotlib.pyplot as plt
from scipy import special 
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
	emc     = -0.583333  # -0.583333
	esc     = -0.416666  # -0.416666
	T_list  = [0.1, 1.0, 10.0]
	colors  = ["darkblue", "gold", "darkred"]
	markers = ["^", "o", "s"]
	alpha   = [1.0, 0.5, 0.5]


	for idx, T in enumerate(T_list):

		chi_ms = 24/T * (ems - 1/2 * emm)
		chi_mc = 24/T * (emc - 1/2 * emm)
		chi_sc = 24/T * (esc) 

		mu_s   = lambda xc: np.log(1 - xc)     + chi_sc * (xc)**2
		mu_c   = lambda xc: np.log(xc)         + chi_sc * (1-xc) ** 2

		xc_bulk  = [1e-6, 0.001, 0.01, 0.1 , 0.2, 0.3, 0.4, 0.4999, 0.5001, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 1-(1e-6)]
		Ns_solv  = []
		Nc_solv  = []
		Nmm_solv = []
		N_solv   = []
		Nmm_mean = 100

		for xc_b in xc_bulk:
			print(f"@ T = {T} and xc in bulk = {xc_b}: mu_s = {mu_s(xc_b)}, mu_c = {mu_c(xc_b)}")
			Nsolv         = lambda Nmm: 770 - 2 * Nmm
			chem_pot      = lambda Nmm, x: -T * Nsolv(Nmm) * (1-x) * mu_s(xc_b) - T * Nsolv(Nmm) * x * mu_c(xc_b)
			entropy_mix   = lambda Nmm, x: -T * Nsolv(Nmm) * ((1-x)*np.log(1-x) + x*np.log(x))
			entropy_chain = lambda Nmm: -T * (Nmm * np.log(Nmm_mean) + (-Nmm_mean) - np.log(special.gamma(Nmm+1)))
			energy        = lambda Nmm, x: (emm * Nmm + ems * Nsolv(Nmm) * (1-x) + emc * Nsolv(Nmm) * x)
			f             = lambda Nmm, x: energy(Nmm, x) + entropy_mix(Nmm, x) + chem_pot(Nmm, x) + entropy_chain(Nmm)

			NMM, XC = np.meshgrid(np.linspace(0,177,1000), np.hstack((0.5, np.linspace(0,1,1000))))
			NMM     = NMM.flatten()
			XC      = XC.flatten()

			free_energy = f(NMM, XC) 
			mask        = np.isnan(free_energy)
			free_energy = free_energy[~mask]
			XC          = XC[~mask]
			idx_min     = np.argmin(free_energy)

			print(f"chem pot = {chem_pot(NMM[idx_min], XC[idx_min])}")
			print(f"entropy_mix  = {entropy_mix(NMM[idx_min], XC[idx_min])}")
			print(f"entropy_chain  = {entropy_chain(NMM[idx_min])}")
			print(f"energy   = {energy(NMM[idx_min], XC[idx_min])}")
			print(f"XC(min)  = {XC[idx_min]}")
			print(f"NMM(min)  = {NMM[idx_min]}")

			Nmm = NMM[idx_min]
			N_solv.append(Nsolv(Nmm))
			Ns_solv.append(Nsolv(Nmm)*(1-XC[idx_min]))
			Nc_solv.append(Nsolv(Nmm)*(XC[idx_min]))
			Nmm_solv.append(Nmm)


		ax[0][0].plot(xc_bulk, Nmm_solv, mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[0][1].plot(xc_bulk, Ns_solv,  mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][0].plot(xc_bulk, Nc_solv,  mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][1].plot(xc_bulk, N_solv,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[0][0].legend(loc="upper right", prop={'size':5})
		ax[0][1].legend(loc="upper right", prop={'size':5})
		ax[1][0].legend(loc="upper right", prop={'size':5})
		ax[1][1].legend(loc="upper right", prop={'size':5})
		ax[1][0].set_xlabel("Bulk number fraction of cosolvent", fontsize=6)
		ax[1][1].set_xlabel("Bulk number fraction of cosolvent", fontsize=6)
		ax[0][0].set_ylabel("Monomer-monomer contacts", fontsize=6) 
		ax[0][1].set_ylabel("Monomer-solvent contacts", fontsize=6)
		ax[1][0].set_ylabel("Monomer-cosolvent contacts", fontsize=6)
		ax[1][1].set_ylabel("Total solvent+cosolvent particles", fontsize=6)
		ax[0][0].set_xlim(0,1)
		ax[0][0].set_ylim(0,180)
		ax[0][1].set_xlim(0,1)
		ax[0][1].set_ylim(0,800)
		ax[1][0].set_xlim(0,1)
		ax[1][0].set_ylim(0,800)
		ax[1][1].set_xlim(0,1)
		ax[1][1].set_ylim(400,800)
	
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
