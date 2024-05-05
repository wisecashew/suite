import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":

	T       = 1
	Nm      = 100
	emm     = -1
	ems     = -1
	emc     = -10
	esc     = -1

	chi_ms = 24/T * (ems - 1/2 * emm)
	chi_mc = 24/T * (emc - 1/2 * emm)
	chi_sc = 24/T * (esc) 

	mu_s   = lambda xc: np.log(xc)     + chi_sc * (1 - xc)**2
	mu_c   = lambda xc: np.log(1 - xc) + chi_sc * xc ** 2

	xc_bulk = [0.001, 0.01, 0.1 , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999]
	Ns_solv = []
	Nc_solv = []
	Nmm_solv  = []


	for xc_b in xc_bulk:
		f = lambda xs, xc: chi_ms * (1 - xs - xc) * xs + chi_mc * (1 - xs - xc) * xc + chi_sc * xs * xc + \
		(1 - xs - xc)/Nm * np.log(1 - xs - xc) + xs * np.log(xs) + xc * np.log(xc) - mu_s(xc_b) * xs/T - mu_c(xc_b) * xc/T

		XS, XC = np.meshgrid(np.linspace(0,1,100), np.linspace(0,1,100))
		XS = XS.flatten()
		XC = XC.flatten()
		mask = (XS + XC < 1)
		XS   = XS[mask]
		XC   = XC[mask]

		free_energy = f(XS, XC) 
		mask        = np.isnan(free_energy)
		free_energy = free_energy[~mask]
		XS          = XS[~mask]
		XC          = XC[~mask]
		idx_min     = np.argmin(free_energy)

		print(f"free energy = {free_energy}", flush=True)
		print(f"min(free energy) = {np.min(free_energy)}", flush=True)
		print(f"minima @ (xs, xc) = ({XS[idx_min]}, {XC[idx_min]})", flush=True)

		Ntot = 1 # Nm/(1-XS[idx_min]-XC[idx_min])
		Ns_solv.append(Ntot * XS[idx_min])
		Nc_solv.append(Ntot * XC[idx_min])
		Nmm_solv.append(1-XS[idx_min]-XC[idx_min])

	fig = plt.figure(figsize=(3,3))
	ax  = plt.axes()

	# ax.plot(xc_bulk, Ns_solv, marker='o', mec='k', ls='--', lw=1, c="coral", label="solvent")
	# ax.plot(xc_bulk, Nc_solv, marker='o', mec='k', ls='--', lw=1, c="steelblue", label="cosolvent")
	ax.plot(xc_bulk, Nmm_solv, marker='o', mec='k', ls='--', lw=1, c="gold", label="frac mm")

	ax.legend(loc="upper right")
	fig.savefig(f"solvationshell.png", dpi=1200, bbox_inches="tight")