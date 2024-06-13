import mpmath as mp
from scipy.stats import beta
import matplotlib.pyplot as plt
import time
import argparse
import itertools
import multiprocessing

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",        dest='o',   type=str,   action='store', help="Name of image.")
parser.add_argument ("--EMM",     dest='emm', type=float, action='store', help="EMM value.")
parser.add_argument ("--EMS",     dest='ems', type=float, action='store', help="EMS value.")
parser.add_argument ("--EMC",     dest='emc', type=float, action='store', help="EMC value.")
parser.add_argument ("--ESC",     dest='esc', type=float, action='store', help="ESC value.")
parser.add_argument ("--T",       dest='T',   nargs='+',  type=float, action='store', help="Enter temperature.")
parser.add_argument ("--color",   dest='c',   nargs='+',  type=str,   action='store', help="color of plot.")
parser.add_argument ("--markers", dest='m',   nargs='+',  type=str,   action='store', help="markers of plot.")
args = parser.parse_args()

def grand_canonical(xc, T):
	
	if xc == 0:
		energy        = lambda Nmm: emm * (Nmm) + ems * Nms_tot(Nmm)
		boltzmann     = lambda Nmm: p_chain(Nmm) * mp.exp(-1/(kb * T) * energy(Nmm))
		# print(f"\tBegin counting averages...")
		ave_Nmm     = 0
		ave_Nms     = 0
		ave_Nmc     = 0
		ave_solvtot = 0
		Z           = 0
		for Nmm in range(31, 208):
			Z += boltzmann(Nmm)
			ave_Nmm     += Nmm * boltzmann(Nmm)
			ave_Nms     += Nms_tot(Nmm)      * boltzmann(Nmm)
			ave_solvtot += Nstot_count(Nmm)  * boltzmann(Nmm)

		# print(f"\tZ = {Z}", flush=True)
		ave_Nmm     = ave_Nmm/Z
		ave_Nms     = ave_Nms/Z
		ave_Nmc     = ave_Nmc/Z
		ave_solvtot = ave_solvtot/Z

	elif xc == 1:
		energy        = lambda Nmm: emm * (Nmm) + emc * Nms_tot(Nmm)
		boltzmann     = lambda Nmm: p_chain(Nmm) * mp.exp(-1/(kb * T) * energy(Nmm))
		Z = 0
		# print(f"\tBegin counting averages...")
		ave_Nmm     = 0
		ave_Nms     = 0
		ave_Nmc     = 0
		ave_solvtot = 0
		for Nmm in range(31, 208):
			Z += boltzmann(Nmm)
			ave_Nmm     += Nmm * boltzmann(Nmm)
			ave_Nmc     += Nms_tot(Nmm)     * boltzmann(Nmm)
			ave_solvtot += Nstot_count(Nmm) * boltzmann(Nmm)

		# print(f"\tZ = {Z}", flush=True)
		ave_Nmm     = ave_Nmm/Z
		ave_Nms     = ave_Nms/Z
		ave_Nmc     = ave_Nmc/Z
		ave_solvtot = ave_solvtot/Z

	else:

		mu_s   = (kb * T) * (mp.log(1 - xc) + chi_sc * (xc) ** 2)   
		mu_c   = (kb * T) * (mp.log(xc)     + chi_sc * (1-xc) ** 2) 

		# print(f"@xc = {xc}: mu_s = {mu_s}, mu_c = {mu_c}")

		omega_mix     = lambda Nmm, Nms: mp.factorial(Nstot_count(Nmm))/(mp.factorial(Nstot_count(Nmm) * Nms/Nms_tot(Nmm)) * mp.factorial(Nstot_count(Nmm) * (1 - Nms/Nms_tot(Nmm))))
		chem_pot      = lambda Nmm, Nms: Nstot_count(Nmm) * Nms/Nms_tot(Nmm) * mu_s + Nstot_count(Nmm) * (1 - Nms/Nms_tot(Nmm)) * mu_c
		energy        = lambda Nmm, Nms: emm * (Nmm) + ems * Nms + emc * (Nms_tot(Nmm) - Nms)
		boltzmann     = lambda Nmm, Nms: mp.exp(-1/(kb * T) * (energy(Nmm, Nms) - chem_pot(Nmm, Nms))) * p_chain(Nmm) * omega_mix(Nmm, Nms) 

		# print(f"\tBegin counting averages...")
		ave_Nmm     = 0
		ave_Nms     = 0
		ave_Nmc     = 0
		ave_solvtot = 0
		Z           = 0
		for Nmm in range(31, 208):
			for Nms in range(0, int(Nms_tot(Nmm))):
				Z += boltzmann(Nmm, Nms)
				ave_Nmm     += Nmm * boltzmann(Nmm, Nms)
				ave_Nms     += Nms * boltzmann(Nmm, Nms)
				ave_Nmc     += (Nms_tot(Nmm) - Nms) * boltzmann(Nmm, Nms)
				ave_solvtot += Nstot_count(Nmm)     * boltzmann(Nmm, Nms)
		
		# print(f"\tZ = {Z}", flush=True)
		ave_Nmm     = ave_Nmm/Z
		ave_Nms     = ave_Nms/Z
		ave_Nmc     = ave_Nmc/Z
		ave_solvtot = ave_solvtot/Z

		print(f"@xc = {xc}: energy = {energy(ave_Nmm, ave_Nms)}, chem pot = {chem_pot(ave_Nmm, ave_Nms)}", flush=True)
	
	return [ave_Nmm, ave_Nms, ave_Nmc, ave_solvtot]


if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(6,4), nrows=2, ncols=2)

	# We change the fontsize of minor ticks label 
	ax[0][0].tick_params(axis='both', which='major', labelsize=4)
	ax[0][1].tick_params(axis='both', which='major', labelsize=4)
	ax[1][0].tick_params(axis='both', which='major', labelsize=4)
	ax[1][1].tick_params(axis='both', which='major', labelsize=4)
	start = time.time()

	Nm      = 32
	z       = 26
	kb      = 1
	emm     = args.emm
	ems     = args.ems
	emc     = args.emc
	esc     = args.esc
	colors  = args.c
	markers = args.m
	alpha   = 1.0
	T_range = args.T

	# parameters for beta distribution
	a     =  6.004256058771041
	b     = 27.319457320842343
	loc   = 31
	scale = 208 - 31

	rv = beta(a, b, loc=loc, scale=scale)
	coeffs = [ 7.45932641e-06, -1.09693434e-02,  5.77531041e+00, -8.28790305e+02]
	Nms_tot       = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain       = lambda Nmm: rv.pdf(Nmm)

	print(f"mean = {rv.mean()}")

	Nstot_count   = lambda Nmm: coeffs[0] * Nms_tot(Nmm) ** 3 + coeffs[1] * Nms_tot(Nmm) ** 2 + coeffs[2] * Nms_tot(Nmm) + coeffs[3]

	print(f"Range of T: {T_range}.", flush=True)
	for idx, T in enumerate(T_range):
		print(f"@T = {T}...", flush=True)
		chi_sc = (z-2)/(kb*T) * (esc) 

		print(f"chi_sc = {chi_sc}")
		xcb        = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # [0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1]
		Nmm_list   = []
		Nms_list   = []
		Nmc_list   = []
		Nsolv_list = []

		pool = multiprocessing.Pool(processes=len(xcb))
		results = pool.starmap(grand_canonical, zip(xcb, itertools.repeat(T)))

		for res in results:
			Nmm_list.append(res[0])
			Nms_list.append(res[1])
			Nmc_list.append(res[2])
			Nsolv_list.append(res[3])

		ax[0][0].plot(xcb, Nmm_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha, clip_on=True)
		ax[0][1].plot(xcb, Nms_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha, clip_on=True)
		ax[1][0].plot(xcb, Nmc_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha, clip_on=True)
		ax[1][1].plot(xcb, Nsolv_list, mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha, clip_on=True)

		
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
	ax[0][0].set_ylim(0,210)
	ax[0][1].set_xlim(0,1)
	ax[0][1].set_ylim(0,800)
	ax[1][0].set_xlim(0,1)
	ax[1][0].set_ylim(0,800)
	ax[1][1].set_xlim(0,1)
	ax[1][1].set_ylim(200,800)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.")

	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
