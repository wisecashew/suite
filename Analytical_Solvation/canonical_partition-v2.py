import mpmath as mp
import numpy as np
from scipy.stats import beta
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
	pv      = 1
	pw      = 0.25
	emm_a   = args.emma
	ems_a   = args.emsa
	emm_n   = args.emmn
	ems_n   = args.emsn
	esc     = 0          # -0.416666
	colors  = [args.c]
	markers = ["^"]
	alpha   = [1.0]

	T_range = [0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.5, 0.7, 0.9, 1.0, 1.25, 1.5, 2.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 50.0, 100.0]
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

	emm_boltz     = lambda T: pv * (pw * mp.exp(-1/(kb * T) * emm_a) * emm_a + (1-pw) * mp.exp(-1/(kb * T) * emm_n) * emm_n)/(pw * mp.exp(-1/(kb * T) * emm_a) + (1-pw) * mp.exp(-1/(kb * T) * emm_n)) + (1-pv) * (emm_n)
	ems_boltz     = lambda T: pv * (pw * mp.exp(-1/(kb * T) * ems_a) * ems_a + (1-pw) * mp.exp(-1/(kb * T) * ems_n) * ems_n)/(pw * mp.exp(-1/(kb * T) * ems_a) + (1-pw) * mp.exp(-1/(kb * T) * ems_n)) + (1-pv) * (ems_n)
	Nms_tot       = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain       = lambda Nmm: rv.pdf(Nmm)
	print(f"Range of T: {T_range}.", flush=True)

	for T in T_range:
		energy        = lambda Nmm: emm_boltz(T) * (Nmm) + ems_boltz(T) * Nms_tot(Nmm)
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
