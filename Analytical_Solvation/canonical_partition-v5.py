import mpmath as mp
import pandas as pd
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-o",      dest='o',     type=str,   action='store', help="Name of image."  )
parser.add_argument ("--EMMA",  dest='emma',  type=float, action='store', help="EMMA value.")
parser.add_argument ("--EMMN",  dest='emmn',  type=float, action='store', help="EMMN value.")
parser.add_argument ("--EMSA",  dest='emsa',  type=float, action='store', help="EMSA value.")
parser.add_argument ("--EMSN",  dest='emsn',  type=float, action='store', help="EMSN value.")
parser.add_argument ("--color", dest='c',     type=str,   action='store', help="color of plot.")
parser.add_argument ("--csv",   dest='csv',   type=str,   action='store', help="name of csv file.")
args = parser.parse_args()

def boltzmann(T):
	energy    = lambda Nmm, Nmm_a, Nms_a: mp.exp(-1/(kb * T) * (Nmm_a * emm_a + (Nmm - Nmm_a) * emm_n + Nms_a * ems_a + (Nms_tot(Nmm) - Nms_a) * ems_n))
	boltzmann = lambda Nmm, Nmm_a, Nms_a: p_chain(Nmm) * comb_mm(Nmm, Nmm_a) * comb_ms(Nmm, Nms_a) * energy(Nmm, Nmm_a, Nms_a)
	print(f"@T = {T}...", flush=True)
	Z = 0
	ave_Nmm       = 0
	ave_Nms       = 0
	for Nmm in range(31, 208):
		print(f"\t@ Nmm = {Nmm}...", flush=True)
		for Nmm_a in range(0, Nmm):
			for Nms_a in range(0, Nms_tot(Nmm)):
				Z        += boltzmann(Nmm, Nmm_a, Nms_a)
				ave_Nmm  += Nmm           * boltzmann(Nmm, Nmm_a, Nms_a)
				ave_Nms  += Nms_tot(Nmm)  * boltzmann(Nmm, Nmm_a, Nms_a)
	return [ave_Nmm/Z, ave_Nms/Z]

if __name__=="__main__":

	fig, ax = plt.subplots(figsize=(8,4), nrows=1, ncols=2)

	# We change the fontsize of minor ticks label 
	ax[0].tick_params(axis='both', which='major', labelsize=6)
	ax[1].tick_params(axis='both', which='major', labelsize=6)
	start = time.time()

	info    = dict()
	Nm      = 32
	z       = 26
	kb      = 1
	pw      = 0.25
	emm_a   = args.emma
	ems_a   = args.emsa
	emm_n   = args.emmn
	ems_n   = args.emsn
	esc     = 0          
	colors  = [args.c]
	markers = ["^", "o"]
	alpha   = [1.0]

	T_range    = [0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
	info["T"]  = T_range
	Nmm_list   = []
	Nms_list   = []

	# parameters for beta distribution
	a     =  6.004256058771041
	b     = 27.319457320842343
	loc   = 31
	scale = 208 - 31

	rv = beta(a, b, loc=loc, scale=scale)
	print(f"p(61) = {rv.pdf(61)}", flush=True)
	print(f"p(80) = {rv.pdf(80)}", flush=True)
	print(f"p(100) = {rv.pdf(100)}", flush=True)
	print(f"p(120) = {rv.pdf(120)}", flush=True)
	print(f"p(150) = {rv.pdf(150)}", flush=True)
	print(f"p(200) = {rv.pdf(200)}", flush=True)
	print(f"p(208) = {rv.pdf(208)}", flush=True)
	exit()

	Nms_tot = lambda Nmm: 770 - 2 * (Nmm - (Nm - 1))
	p_chain = lambda Nmm: rv.pdf(Nmm)
	comb_mm = lambda Nmm, Nmm_a: mp.factorial(Nmm)          / (mp.factorial(Nmm_a) * mp.factorial(Nmm-Nmm_a))          * mp.power(pw, Nmm_a) * mp.power(1-pw, Nmm-Nmm_a)
	comb_ms = lambda Nmm, Nms_a: mp.factorial(Nms_tot(Nmm)) / (mp.factorial(Nms_a) * mp.factorial(Nms_tot(Nmm)-Nms_a)) * mp.power(pw, Nms_a) * mp.power(1-pw, Nms_tot(Nmm)-Nms_a)

	pool    = multiprocessing.Pool(processes=len(T_range))
	results = pool.starmap(boltzmann, zip(T_range))

	for res in results:
		Nmm_list.append(res[0])
		Nms_list.append(res[1])

	info["Nmm"] = Nmm_list
	info["Nms"] = Nms_list

	df = pd.DataFrame(info)
	df.to_csv(args.csv)

	ax[0].set_xscale('log')
	ax[1].set_xscale('log')
	ax[0].plot(T_range, Nmm_list, mec='k', ls='--', lw=1, c=colors[0], marker=markers[0], markersize=8, alpha=alpha[0])
	ax[1].plot(T_range, Nms_list, mec='k', ls='--', lw=1, c=colors[0], marker=markers[1], markersize=8, alpha=alpha[0])
	ax[0].set_xlabel("Temperature (T)", fontsize=8)
	ax[1].set_xlabel("Temperature (T)", fontsize=8)
	ax[0].set_ylabel("MM contacts", fontsize=8)
	ax[1].set_ylabel("MS contacts", fontsize=8)
	ax[0].set_xlim(np.min(T_range), np.max(T_range))
	ax[0].set_ylim(0,210)
	ax[1].set_xlim(np.min(T_range), np.max(T_range))
	ax[1].set_ylim(400,800)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")
