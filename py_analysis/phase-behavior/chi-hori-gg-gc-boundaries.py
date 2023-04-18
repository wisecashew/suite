#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator
import scipy.optimize as so
import copy
from scipy.stats import linregress



if __name__=="__main__":


    ems_list = [-1.4,-1.5]
    # elow    = -50
    # ehigh   = -1.4
    g  = 0.5
    
    zmm  = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float64) + (1-g)*np.exp ((-1/T * emmn), dtype = np.float64)
    zms  = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float64) + (1-g)*np.exp ((-1/T * emsn), dtype = np.float64)
    fmma = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float64) / zmm(emma, emmn, T)
    fmsa = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float64) / zms(emsa, emsn, T)

    def chi (emma, emmn, emsa, emsn, pv, T):
        t1 = pv*(fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn) + (1-pv)*emmn
        return 1/T * (t1 - 0.5 * t2)   


    def dfmma (emma, emmn, T):

        df = 1/(T**2) * ( np.exp((emma+emmn)/T, dtype=np.float128) * (emma - emmn) * g * (1-g) ) / (g * np.exp(emmn/T, dtype=np.float128) + (1-g)*np.exp(emma/T, dtype=np.float128))**2

        return df 

    def dfmsa (emsa, emsn, T):

        df = 1/(T**2) * ( np.exp((emsa+emsn)/T, dtype=np.float128) * (emsa - emsn) * g * (1-g) ) / (g * np.exp (emsn/T, dtype=np.float128) + (1-g) * np.exp(emsa/T, dtype=np.float128) )**2

        return df 

    def dchi (emma, emmn, emsa, emsn, pv, T):

        term1 = -1/T * chi (emma, emmn, emsa, emsn, pv, T) 
        term2 =  1/T * ( dfmsa (emsa, emsn, T) * (emsa-emsn) - 1/2 * dfmma (emma, emmn, T) *(emma-emmn) )

        return term1 + term2


    def check_temp_interval (E_list, pv):

        sgn = np.sign( chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 0.01) ) * np.sign(chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 100))

        if sgn == 1:
            return 0

        else:
            return 1

    def bisection (E_list, pv, interval, itermax):

        for i in range(itermax):

            E_list_l    = copy.copy(E_list) 
            E_list_l[1] = interval[0]
            E_list_r    = copy.copy(E_list)
            E_list_r[1] = interval[1]

            res_l = check_temp_interval (E_list_l, pv) 
            res_r = check_temp_interval (E_list_r, pv)

            if res_l + res_r != 1:
                print ("Bad bounds!")
                exit()
                return None

            else:
                E_list_m = copy.copy(E_list)
                E_list_m[1] = np.mean(interval)
                res_m = check_temp_interval (E_list_m, pv)

            if res_l + res_m == 1:
                interval[1] = np.mean(interval)

            else:
                interval[0] = np.mean(interval)

            chi_l = lambda T: chi (E_list[0], interval[1], E_list[2], E_list[3], pv, T)
            sol = so.brentq (chi_l, 0.01, 100)

            if np.abs(100-sol) < 1e-4:
                return interval[1]

            else:
                E_list[1] = interval[1]
        
        return interval[1]
    



    norm = matplotlib.colors.TwoSlopeNorm (vmin=np.min(ems_list), vcenter=np.mean(ems_list), vmax=np.max(ems_list))

    E_ms_a = -1   ;
    E_ms_n = -0   ;
    pv     = 1.0  ;
    
    y = np.arange (-3.5, -2, 0.01)
    x = []

    interval = [-2.5, 2.5]
    E_list   = [0, 0, E_ms_a, E_ms_n]
    for E_mm_a in np.arange (-3.5, -2, 0.01):
        E_list [0] = E_mm_a
        interval = [-2.5, 2.5]
        E_mm_n = bisection (E_list, pv, interval, 1000)
        x.append (E_mm_n)
    
    print (f"E_mm_n = {x}")
    slope, intercept, r_val, p_val, std_err = linregress (x, y)
    r2 = r_val ** 2

    print (f"Slope: {slope}")
    print (f"intercept: {intercept}")
    print (f"r2: {r2}")


    # for E_mm_n in [0.537]: # np.linspace (-2.5, 2.5, 1000):
    #    chi_l = lambda T: chi (E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, T)
    #    sol = so.brentq (chi_l, 0.01, 100)
    #     print (f"sol = {sol}")
    #     print (f"chi = {chi_l(sol)}")
    #     print (f"dchi = {dchi(E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, sol)}")

    
    

    
    
