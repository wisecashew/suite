import mpmath as mp
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",      dest='o',     type=str, action='store', help="Name of image."  )
parser.add_argument ("--EMMA",  dest='emma',  type=float, action='store', help="EMMA value.")
parser.add_argument ("--EMMN",  dest='emmn',  type=float, action='store', help="EMMN value.")
parser.add_argument ("--EMSA",  dest='emsa',  type=float, action='store', help="EMSA value.")
parser.add_argument ("--EMSN",  dest='emsn',  type=float, action='store', help="EMSN value.")
parser.add_argument ("--color", dest='c',     type=str, action='store', help="color of plot.")
# parser.add_argument ("--nproc", dest='nproc', type=int, action='store', help="number of processors.")
args = parser.parse_args()

def boltzmann(T):
	energy    = lambda Nmm, Nmm_a, Nms_a: mp.exp(-1/(kb * T) * (Nmm_a * emm_a + (Nmm - Nmm_a) * emm_n + Nms_a * ems_a + (Nms_tot(Nmm) - Nms_a) * ems_n))
	boltzmann = lambda Nmm, Nmm_a, Nms_a: p_chain(Nmm) * comb_mm(Nmm, Nmm_a) * comb_ms(Nmm, Nms_a) * energy(Nmm, Nmm_a, Nms_a)
	print(f"@T = {T}...", flush=True)
	Z = 0
	ave_Nmm       = 0
	ave_Nms       = 0
	ave_Nmm_store = []
	ave_Nms_store = []
	for Nmm in range(31, 208):
		print(f"\t@Nmm = {Nmm}...", flush=True)
		for Nmm_a in range(0, Nmm):
			for Nms_a in range(0, Nms_tot(Nmm)):
				Z += boltzmann(Nmm, Nmm_a, Nms_a)
				ave_Nmm     += Nmm           * boltzmann(Nmm, Nmm_a, Nms_a)
				ave_Nms     += Nms_tot(Nmm)  * boltzmann(Nmm, Nmm_a, Nms_a)
				ave_Nmm_store.append(ave_Nmm)
				ave_Nms_store.append(ave_Nms)
	ave_Nmm_store = np.array(ave_Nmm_store)
	ave_Nms_store = np.array(ave_Nms_store)
	ave_Nmm = ave_Nmm/Z
	ave_Nms = ave_Nms/Z
	print(f"Max Nmm term = {np.max(ave_Nmm_store)/Z}, avg = {ave_Nmm}, T = {T}", flush=True)
	print(f"Max Nms term = {np.max(ave_Nmm_store)/Z}, avg = {ave_Nms}, T = {T}", flush=True)


	# print(f"\tZ = {Z}", flush=True)
	return [ave_Nmm, ave_Nms]

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

	# parameters for beta distribution
	a     = 6.004256058771041  # 8.382763329911416
	b     = 27.319457320842343 # 1944745.1084205871
	loc   = 31       # 28.403244859804218
	scale = 208 - 31 # 7751700.2611331325

	rv = beta(a, b, loc=loc, scale=scale)

	Nms_tot = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain = lambda Nmm: rv.pdf(Nmm)
	comb_mm = lambda Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(pw, Nmm_a) * mp.power(1-pw, Nmm-Nmm_a)
	comb_ms = lambda Nmm, Nms_a: mp.factorial(Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(Nms_tot(Nmm)-Nms_a)) * mp.power(pw, Nms_a) * mp.power(1-pw, Nms_tot(Nmm)-Nms_a)

	print(f"Range of T: {T_range}.", flush=True)
	pool    = multiprocessing.Pool(processes=len(T_range))
	results = pool.starmap(boltzmann, zip(T_range))

	for res in results:
		Nmm_list.append(res[0])
		Nms_list.append(res[1])


	print(f"Nmm_list = {Nmm_list}", flush=True)
	print(f"Nms_list = {Nms_list}", flush=True)

	ax[0].set_xscale('log')
	ax[1].set_xscale('log')
	ax[0].plot(T_range, Nmm_list, mec='k', ls='--', lw=1, c=colors[0], marker=markers[0], alpha=alpha[0])
	ax[1].plot(T_range, Nms_list, mec='k', ls='--', lw=1, c=colors[0], marker=markers[0], alpha=alpha[0])
	ax[1].set_xlabel("Temperature (T)", fontsize=8)
	ax[0].set_ylabel("MM contacts", fontsize=8)
	ax[1].set_ylabel("MS contacts", fontsize=8)
	ax[0].set_xlim(np.min(T_range), np.max(T_range))
	ax[0].set_ylim(0,210)
	ax[1].set_xlim(np.min(T_range), np.max(T_range))
	ax[1].set_ylim(400,800)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.")

	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
