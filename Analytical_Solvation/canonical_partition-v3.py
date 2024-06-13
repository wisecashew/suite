import mpmath as mp
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",      dest='o', type=str, action='store', help="Name of image."  )
parser.add_argument ("--EMMA",  dest='emma', type=float, action='store', help="EMMA value.")
parser.add_argument ("--EMMN",  dest='emmn', type=float, action='store', help="EMMN value.")
parser.add_argument ("--EMSA",  dest='emsa', type=float, action='store', help="EMSA value.")
parser.add_argument ("--EMSN",  dest='emsn', type=float, action='store', help="EMSN value.")
parser.add_argument ("--color", dest='c',   type=str, action='store', help="color of plot.")
args = parser.parse_args()


if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(4,4), nrows=2, ncols=1)

	# We change the fontsize of minor ticks label 
	ax[0].tick_params(axis='both', which='major', labelsize=6)
	ax[1].tick_params(axis='both', which='major', labelsize=6)
	start = time.time()

	Nm      = 32
	z       = 26
	kb      = 1
	pw      = 0.25
	emm_a   = args.emma
	ems_a   = args.emsa
	emm_n   = args.emmn
	ems_n   = args.emsn
	esc     = 0          # -0.416666
	colors  = [args.c]
	markers = ["^"]
	alpha   = [1.0]

	T_range    = [0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
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

	Nms_tot = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain = lambda Nmm: rv.pdf(Nmm)
	comb_mm = lambda Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(pw, Nmm_a) * mp.power(1-pw, Nmm-Nmm_a)
	comb_ms = lambda Nmm, Nms_a: mp.factorial(Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(Nms_tot(Nmm)-Nms_a)) * mp.power(pw, Nms_a) * mp.power(1-pw, Nms_tot(Nmm)-Nms_a)

	print(f"Range of T: {T_range}.", flush=True)

	for T in T_range:
		energy    = lambda Nmm, Nmm_a, Nms_a: mp.exp(-1/(kb * T) * (Nmm_a * emm_a + (Nmm - Nmm_a) * emm_n + Nms_a * ems_a + (Nms_tot(Nmm) - Nms_a) * ems_n))
		boltzmann = lambda Nmm, Nmm_a, Nms_a: p_chain(Nmm) * comb_mm(Nmm, Nmm_a) * comb_ms(Nmm, Nms_a) * energy(Nmm, Nmm_a, Nms_a)
		print(f"@T = {T}...", flush=True)
		Z = 0
		for Nmm in range(31, 208):
			print(f"\t@Nmm = {Nmm}", flush=True)
			for Nmm_a in range(0, Nmm):
				for Nms_a in range(0, Nms_tot(Nmm)):
					Z += boltzmann(Nmm, Nmm_a, Nms_a)

		print(f"\tZ = {Z}", flush=True)
		print(f"\tBegin counting averages...")
		ave_Nmm     = 0
		ave_Nms     = 0

		for Nmm in range(31, 208):
			for Nmm_a in range(0, Nmm):
				for Nms_a in range(0, Nms_tot(Nmm)):
					ave_Nmm     += Nmm * boltzmann(Nmm, Nmm_a, Nms_a)
					ave_Nms     += Nms_tot(Nmm)  * boltzmann(Nmm, Nmm_a, Nms_a)/Z

		Nmm_list.append(ave_Nmm)
		Nms_list.append(ave_Nms)
		print(f"\tave_Nmm = {ave_Nmm}, ave_Nms = {ave_Nms}.", flush=True)
		
		print(f"Done with averaging. Plotting...")
		idx = 0

	ax[0].plot(T_range, Nmm_list, mec='k', ls='--', lw=1, c=colors[idx], marker=markers[idx], alpha=alpha[idx], clip_on=False)
	ax[1].plot(T_range, Nms_list, mec='k', ls='--', lw=1, c=colors[idx], marker=markers[idx], alpha=alpha[idx], clip_on=False)
	ax[1].set_xlabel("Temperature (T)", fontsize=8)
	ax[0].set_ylabel("MM contacts", fontsize=8)
	ax[1].set_ylabel("MS contacts", fontsize=8)
	ax[0].set_xlim(np.min(T_range), np.max(T_range))
	ax[0].set_ylim(0,210)
	ax[1].set_xlim(np.min(T_range), np.max(T_range))
	ax[1].set_ylim(400,800)
	ax[0].set_xscale('log')
	ax[1].set_xscale('log')

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.")

	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
