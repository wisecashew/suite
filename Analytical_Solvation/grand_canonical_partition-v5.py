import mpmath as mp
import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",      dest='o',   type=str,   action='store', help="Name of image.")
parser.add_argument ("--EMM",   dest='emm', type=float, action='store', help="EMM value.")
parser.add_argument ("--EMS",   dest='ems', type=float, action='store', help="EMS value.")
parser.add_argument ("--EMC",   dest='emc', type=float, action='store', help="EMC value.")
parser.add_argument ("--ESC",   dest='esc', type=float, action='store', help="ESC value.")
parser.add_argument ("--T",     dest='T',   type=float, action='store', help="Enter temperature.")
parser.add_argument ("--color", dest='c',   type=str,   action='store', help="color of plot.")
args = parser.parse_args()

if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(6,4), nrows=2, ncols=2)

	# We change the fontsize of minor ticks label 
	ax[0][0].tick_params(axis='both', which='major', labelsize=4)
	ax[0][1].tick_params(axis='both', which='major', labelsize=4)
	ax[1][0].tick_params(axis='both', which='major', labelsize=4)
	ax[1][1].tick_params(axis='both', which='major', labelsize=4)
	start = time.time()

	Nm       = 32
	z        = 26
	kb       = 1
	emm      = args.emm
	ems      = args.ems
	emc      = args.emc
	esc      = args.esc
	colors   = [args.c]
	markers  = ["^"]
	alpha    = [1.0]
	T_range  = [args.T]
	L        = 34
	Tsolvent = L ** 3 - Nm

	print(f"Range of T: {T_range}.", flush=True)
	for T in T_range:

		chi_sc = (z-2)/(kb*T) * (esc) 

		xcb        = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
		Nmm_list   = []
		Nms_list   = []
		Nmc_list   = []
		Nsolv_list = []

		# parameters for beta distribution
		a     = 6.004256058771041  # 8.382763329911416
		b     = 27.319457320842343 # 1944745.1084205871
		loc   = 31       # 28.403244859804218
		scale = 208 - 31 # 7751700.2611331325

		rv = beta(a, b, loc=loc, scale=scale)
		coeffs = [ 7.45932641e-06, -1.09693434e-02,  5.77531041e+00, -8.28790305e+02]
		print(f"mean = {rv.mean()}")
		print(f"@T = {T}...", flush=True)

		# this is the total number of unique particles in the solvation shell 
		Nstot_count   = lambda Nmm: coeffs[0] * Nms_tot(Nmm) ** 3 + coeffs[1] * Nms_tot(Nmm) ** 2 + coeffs[2] * Nms_tot(Nmm) + coeffs[3]

		for xc_bulk in xcb:
			print(f"xc = {xc_bulk}...", flush=True)
			Nms_tot       = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
			p_chain       = lambda Nmm: rv.pdf(Nmm)
			Ns            = int(Tsolvent * (1 - xc_bulk))
			Nc            = Tsolvent - Ns - Nm

			if xc_bulk == 0:
				energy        = lambda Nmm: emm * (Nmm) + ems * Nms_tot(Nmm)
				boltzmann     = lambda Nmm: p_chain(Nmm) * mp.exp(-1/(kb * T) * energy(Nmm))
				print(f"\tBegin counting averages...")
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

				print(f"\tZ = {Z}", flush=True)
				ave_Nmm     = ave_Nmm/Z
				ave_Nms     = ave_Nms/Z
				ave_Nmc     = ave_Nmc/Z
				ave_solvtot = ave_solvtot/Z

			elif xc_bulk == 1:
				energy        = lambda Nmm: emm * (Nmm) + emc * Nms_tot(Nmm)
				boltzmann     = lambda Nmm: p_chain(Nmm) * mp.exp(-1/(kb * T) * energy(Nmm))
				Z = 0
				print(f"\tBegin counting averages...")
				ave_Nmm     = 0
				ave_Nms     = 0
				ave_Nmc     = 0
				ave_solvtot = 0
				for Nmm in range(31, 208):
					Z += boltzmann(Nmm)
					ave_Nmm     += Nmm * boltzmann(Nmm)
					ave_Nmc     += Nms_tot(Nmm)     * boltzmann(Nmm)
					ave_solvtot += Nstot_count(Nmm) * boltzmann(Nmm)

				print(f"\tZ = {Z}", flush=True)
				ave_Nmm     = ave_Nmm/Z
				ave_Nms     = ave_Nms/Z
				ave_Nmc     = ave_Nmc/Z
				ave_solvtot = ave_solvtot/Z

			else:
				N_min      = np.min([Ns, Nc])
				energy     = lambda Nmm, Nms, Nsc: emm * (Nmm) + ems * Nms + emc * (Nms_tot(Nmm) - Nms) + esc * Nsc
				boltzmann  = lambda Nmm, Nms, Nsc: p_chain(Nmm) * mp.exp(-1/(kb * T) * (energy(Nmm, Nms, Nsc)))

				print(f"\tBegin counting averages...")
				ave_Nmm     = 0
				ave_Nms     = 0
				ave_Nmc     = 0
				ave_solvtot = 0
				Z           = 0
				for Nmm in range(31, 208):
					print(f"Nmm = {Nmm}", flush=True)
					for Nms in range(0, int(Nms_tot(Nmm))):
						print(f"\tNms = {Nms}", flush=True)
						Nsc_max = np.min([936360-Nms+z*int(39272*(1-xc_bulk)), 1020240+2*Nmm+Nms-26*int(39272*(1-xc_bulk))])
						Ncube   = int(np.floor(N_min ** (1/3)))
						Nsc_min = 54 * (Ncube-2)**2 + 12 * 15 * (Ncube-2) + 19*8
						print(f"\t\tNsc range: [{Nsc_min}, {Nsc_max}]")
						for Nsc in range(Nsc_min, Nsc_max):
							Z           += boltzmann(Nmm, Nms, Nsc)
							ave_Nmm     += Nmm * boltzmann(Nmm, Nms, Nsc)
							ave_Nms     += Nms * boltzmann(Nmm, Nms, Nsc)
							ave_Nmc     += (Nms_tot(Nmm) - Nms) * boltzmann(Nmm, Nms, Nsc)
							ave_solvtot += Nstot_count(Nmm)     * boltzmann(Nmm, Nms, Nsc)
				
				print(f"\tZ = {Z}", flush=True)
				ave_Nmm     = ave_Nmm/Z
				ave_Nms     = ave_Nms/Z
				ave_Nmc     = ave_Nmc/Z
				ave_solvtot = ave_solvtot/Z

			Nmm_list.append  (ave_Nmm)
			Nms_list.append  (ave_Nms)
			Nmc_list.append  (ave_Nmc)
			Nsolv_list.append(ave_solvtot)
			print(f"\tave_Nmm = {ave_Nmm}, ave_Ns = {ave_Nms}, ave_Nc = {ave_Nmc}, ave_Nsolv = {ave_solvtot}.", flush=True)
		
		print(f"Done with averaging. Plotting...")
		idx = 0
		ax[0][0].plot(xcb, Nmm_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[0][1].plot(xcb, Nms_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][0].plot(xcb, Nmc_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][1].plot(xcb, Nsolv_list, mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)

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
