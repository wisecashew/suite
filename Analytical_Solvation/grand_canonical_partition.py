import mpmath as mp
from scipy.stats import beta
import matplotlib.pyplot as plt
import time

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
	emm     = -1
	ems     = -0.541666
	emc     = -0.916666  # -0.583333
	esc     = 0          # -0.416666
	colors  = ["darkblue", "gold", "darkred"]
	markers = ["^", "o", "s"]
	alpha   = [1.0, 0.5, 0.5]

	T_range = [10.0]
	print(f"Range of T: {T_range}.", flush=True)
	for T in T_range:
		
		chi_ms = (z-2)/(kb*T) * (ems - 1/2 * emm)
		chi_mc = (z-2)/(kb*T) * (emc - 1/2 * emm)
		chi_sc = (z-2)/(kb*T) * (esc) 

		xcb = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
		Nmm_list   = []
		Ns_list    = []
		Nc_list    = []
		Nsolv_list = []
		# parameters for beta distribution
		a=8.382763329911416
		b=1944745.1084205871
		loc=28.403244859804218
		scale=7751700.2611331325

		rv = beta(a, b, loc=loc, scale=scale)
		coeffs = [ 7.45932641e-06, -1.09693434e-02,  5.77531041e+00, -8.28790305e+02]
		print(f"mean = {rv.mean()}")
		print(f"@T = {T}...", flush=True)
		for xc_bulk in xcb:

			mu_s   = (kb * T) * (mp.log(1 - xc_bulk) + chi_sc * (xc_bulk) ** 2)   if xc_bulk > 0 else 0
			mu_c   = (kb * T) * (mp.log(xc_bulk)     + chi_sc * (1-xc_bulk) ** 2) if xc_bulk > 0 else 0

			print(f"\t@xc in bulk = {xc_bulk}...", end=' ', flush=True) 
			Nms_tot       = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
			p_chain       = lambda Nmm: rv.pdf(Nmm)
			Nstot_count   = lambda Nmm: coeffs[0] * Nms_tot(Nmm) ** 3 + coeffs[1] * Nms_tot(Nmm) ** 2 + coeffs[2] * Nms_tot(Nmm) + coeffs[3]
			omega_mix     = lambda Nmm, Ns: mp.factorial(Nstot_count(Nmm))/(mp.factorial(Ns) * mp.factorial(Nstot_count(Nmm)-Ns))
			chem_pot      = lambda Nmm, Ns: Ns * mu_s + (Nstot_count(Nmm) - Ns) * mu_c
			energy        = lambda Nmm, Ns: emm * (Nmm) + ems * Nms_tot(Nmm) * Ns/Nstot_count(Nmm) + emc * Nms_tot(Nmm) * (1 - Ns/Nstot_count(Nmm))
			boltzmann     = lambda Nmm, Ns: p_chain(Nmm) * omega_mix(Nmm, Ns) * mp.exp(-1/(kb * T) * (energy(Nmm, Ns) - chem_pot(Nmm, Ns)))

			Z = 0
			for Nmm in range(31, 208):
				for Ns in range(0, int(Nstot_count(Nmm))):
					Z += boltzmann(Nmm, Ns)

			print(f"\tZ = {Z}", flush=True)
			print(f"\tBegin counting averages...")
			ave_Nmm     = 0
			ave_Ns      = 0
			ave_Nc      = 0
			ave_solvtot = 0

			for Nmm in range(31, 208):
				for Ns in range(0, int(Nstot_count(Nmm))):
					ave_Nmm     += Nmm * boltzmann(Nmm, Ns)/Z
					ave_Ns      += Ns  * boltzmann(Nmm, Ns)/Z
					ave_Nc      += (Nstot_count(Nmm) - Ns) * boltzmann(Nmm, Ns)/Z
					ave_solvtot += Nstot_count(Nmm) * boltzmann(Nmm, Ns)/Z

			Nmm_list.append(ave_Nmm)
			Ns_list.append(ave_Ns)
			Nc_list.append(ave_Nc)
			Nsolv_list.append(ave_solvtot)
			print(f"\tave_Nmm = {ave_Nmm}, ave_Ns = {ave_Ns}, ave_Nc = {ave_Nc}, ave_Nsolv = {ave_solvtot}.", flush=True)
		
		print(f"Done with averaging. Plotting...")
		idx = 0
		ax[0][0].plot(xcb, Nmm_list,     mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[0][1].plot(xcb, Ns_list,      mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][0].plot(xcb, Nc_list,      mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
		ax[1][1].plot(xcb, Nsolv_list,   mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)

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
	ax[1][1].set_ylim(300,800)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.")

	fig.savefig("averages", dpi=1200, bbox_inches="tight")
