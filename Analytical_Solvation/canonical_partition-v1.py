import mpmath as mp
from scipy.stats import beta
import matplotlib.pyplot as plt
import time

if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(4,4), nrows=2, ncols=1)

	# We change the fontsize of minor ticks label 
	ax[0].tick_params(axis='both', which='major', labelsize=4)
	ax[1].tick_params(axis='both', which='major', labelsize=4)
	start = time.time()

	Nm      = 32
	z       = 26
	kb      = 1
	emm     = -1.5
	ems     = -0.541666
	emc     = -0.916666  # -0.583333
	esc     = 0          # -0.416666
	colors  = ["darkblue", "gold", "darkred"]
	markers = ["^", "o", "s"]
	alpha   = [1.0, 0.5, 0.5]

	T_range = [0.1, 0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 50.0, 100.0]
	Nmm_list   = []
	Nms_list   = []
	Nc_list    = []
	Nsolv_list = []

	# parameters for beta distribution
	a     = 8.382763329911416
	b     = 1944745.1084205871
	loc   = 28.403244859804218
	scale = 7751700.2611331325

	rv = beta(a, b, loc=loc, scale=scale)

	Nms_tot       = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain       = lambda Nmm: rv.pdf(Nmm)
	energy        = lambda Nmm: emm * (Nmm) + ems * Nms_tot(Nmm)
	
	print(f"Range of T: {T_range}.", flush=True)
	for T in T_range:

		boltzmann     = lambda Nmm: p_chain(Nmm) * mp.exp(-1/(kb * T) * (energy(Nmm)))	
		print(f"@T = {T}...", flush=True)
		Z = 0
		for Nmm in range(31, 208):
			Z += boltzmann(Nmm)

		print(f"\tZ = {Z}", flush=True)
		print(f"\tBegin counting averages...")
		ave_Nmm     = 0
		ave_Nms     = 0

		for Nmm in range(31, 208):
			ave_Nmm     += Nmm * boltzmann(Nmm)/Z
			ave_Nms     += Nms_tot(Nmm)  * boltzmann(Nmm)/Z

		Nmm_list.append(ave_Nmm)
		Nms_list.append(ave_Nms)
		print(f"\tave_Nmm = {ave_Nmm}, ave_Nms = {ave_Nms}.", flush=True)
		
		print(f"Done with averaging. Plotting...")
		idx = 0

	ax[0].plot(T_range, Nmm_list, mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)
	ax[1].plot(T_range, Nms_list, mec='k', ls='--', lw=1, c=colors[idx], label=f"T={T}", marker=markers[idx], alpha=alpha[idx], clip_on=False)

	ax[0].legend(loc="upper right", prop={'size':5})
	ax[1].legend(loc="upper right", prop={'size':5})
	ax[0].set_ylabel("Monomer-monomer contacts", fontsize=6) 
	ax[1].set_ylabel("Monomer-solvent contacts", fontsize=6)
	ax[0].set_xlim(0.1,100)
	ax[0].set_ylim(0,210)
	ax[1].set_xlim(0.1,100)
	ax[1].set_ylim(400,800)
	ax[0].set_xscale('log')
	ax[1].set_xscale('log')

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.")

	fig.savefig("single-solvent", dpi=1200, bbox_inches="tight")
