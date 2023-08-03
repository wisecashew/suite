import numpy as np
import matplotlib.pyplot as plt 

if __name__=="__main__":

    chi_eff  = lambda beta, delta_ms, delta_mm, e_ms_n, e_mm_n, z, r: beta*(z-2)*(np.exp (-beta*delta_ms) * delta_ms/(np.exp(-beta*delta_ms)+r) - 1/2 * np.exp(-beta*delta_mm)*delta_mm/(np.exp(-beta*delta_mm)+r) + e_ms_n - 1/2 * e_mm_n) - 0.692

    def beta_sol (delta_ms, delta_mm, e_ms_n, e_mm_n, z, r):

        d1_mm = 1 + np.exp (delta_mm)*r
        d1_ms = 1 + np.exp (delta_ms)*r
        t1    = np.exp(delta_mm)*(delta_mm)**2/d1_mm**2 - 1/2*np.exp(delta_mm)*z*(delta_mm)**2/d1_mm**2 + 0.5*np.exp(delta_mm)*delta_mm**3/()

