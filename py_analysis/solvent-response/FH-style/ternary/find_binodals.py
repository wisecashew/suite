import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import argparse

if __name__=="__main__":

	# binodal plotter
    df = pd.read_csv ("comparing_mu.txt", sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])

    df = df.loc[df["dmu"]<0.01]

    # print (df)

    va = 1
    vb = 32
    vc = 1

    chi_ab = -1
    chi_ac = 0
    chi_bc = -10

    phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
    phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
    phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    for idx in range (1, 300, 10):

	    def mu_equations (phi):
	    	eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
	    	eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
	    	eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

	    	return [eq1, eq2, eq3]


	    root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]])
	    print (f"phi_b[{idx}] = {phi_b[idx]}")
	    print (f"roots = {root}")
	    print (f"delta mu = {mu_equations(root)}")