#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator
import scipy.optimize as so
import copy
from scipy.stats import linregress
import time
import pandas as pd
from sklearn.cluster import KMeans

#######################################################


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
    T   = np.logspace(-2, 2, 100)
    sgn = np.sign( chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 0.01) ) * np.sign(np.max (chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, T) ) )

    if sgn == 1: # if same sign
        return 0

    else: # if opposite sign
        return 1


def range_finder (E_mm_a_range, E_mm_n_range, E_ms_list, pv):

    no_overlap   = []
    with_overlap = []


    for E_mm_a in E_mm_a_range:
        E_list = [E_mm_a, E_mm_n_range[0 ], E_ms_list[0], E_ms_list[1]]
        c1 = check_temp_interval (E_list, pv)
        E_list = [E_mm_a, E_mm_n_range[-1], E_ms_list[0], E_ms_list[1]]
        c2 = check_temp_interval (E_list, pv)

        if c1 + c2 == 1:
            no_overlap.append   (E_mm_a)
        else:
            with_overlap.append (E_mm_a)

    return no_overlap, with_overlap


def bisection (E_list, pv, interval, itermax):

    for i in range(itermax):

        E_list_l    = copy.copy(E_list) 
        E_list_l[1] = interval[0]
        E_list_r    = copy.copy(E_list)
        E_list_r[1] = interval[1]

        res_l = check_temp_interval (E_list_l, pv) 
        res_r = check_temp_interval (E_list_r, pv)

        if res_l + res_r != 1:
            print (f"Low T sign = {np.sign( chi (E_list_l[0], E_list_l[1], E_list_l[2], E_list_l[3], pv, 0.01) )}")
            print (f"High T sign = {np.sign( np.max ( chi (E_list_l[0], E_list_l[1], E_list_l[2], E_list_l[3], pv,  np.logspace(-2,2,100) ) ) ) }")
            print (f"E_list_l = {E_list_l}")
            print (f"E_list_r = {E_list_r}")
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

        sol   = np.max (chi (E_list[0], interval[0], E_list[2], E_list[3], pv, np.logspace(-2,2,100) ) ) # so.brentq (chi_l, 0.01, 100)

        if np.abs(sol) < 1e-4:
            return interval[0]

        else:
            E_list[1] = interval[0]
    
    return interval[0]


def line_search (E_list, E_mm_n_range, pv):
    
    line_bools = []
    for E_mm_n in E_mm_n_range:
        E_list[1] = E_mm_n
        line_bools.append ( check_temp_interval (E_list, pv) )

    line_bools = np.array (line_bools)
    line_bools = line_bools==0
    if len(line_bools[line_bools]) == 0:
        return [False]
    else:
        E_left  = E_mm_n_range[line_bools][0 ]
        E_right = E_mm_n_range[line_bools][-1]
        return [True, (E_left, E_right)]



if __name__=="__main__":

    start = time.time()
    E_ms_a = -1   ;
    E_ms_n = -0   ;
    pv     = 1.0  ;
    
    E_mm_a_range = np.arange (-1.99, 2.5, 0.001)
    E_mm_n_range = np.arange (-1.99, 2.5, 0.001)
    E_ms_list = [E_ms_a, E_ms_n]

    E_mm_a_no_overlap, E_mm_a_with_overlap = range_finder (E_mm_a_range, E_mm_n_range, E_ms_list, pv)

    E_mm_n_list = []

    # for points with no overlap, I can use bisection 
    for E_mm_a in E_mm_a_no_overlap:
        E_list = [E_mm_a, 0, E_ms_a, E_ms_n]
        E_mm_n_interval = [-1.9, 2.5]
        E_mm_n_list.append ( bisection (E_list, pv, E_mm_n_interval, 100) )


    data = {"E_mm_a": E_mm_a_no_overlap, "E_mm_n": E_mm_n_list}
    df = pd.DataFrame.from_dict (data)
    df.to_csv ("boundaries_tophalf.csv", index=False)


    left_arm  = []
    right_arm = []
    E_mm_a_wo_plot = []
    for E_mm_a in E_mm_a_with_overlap:
        E_list = [E_mm_a, 0, E_ms_a, E_ms_n]
        results = line_search (E_list, E_mm_n_range, pv)

        if results[0]:
            left_arm.append  (results[1][0])
            right_arm.append (results[1][1])
            E_mm_a_wo_plot.append (E_mm_a)
        else:
            continue


    data = {"E_mm_a": E_mm_a_wo_plot, "E_mm_n": right_arm}
    df   = pd.DataFrame.from_dict (data)
    df.to_csv ("boundaries_rightarm.csv", index=False)


    plt.figure(1)
    plt.plot (E_mm_n_list, E_mm_a_no_overlap, color='steelblue')
    plt.plot (left_arm , E_mm_a_wo_plot, color='darkred')
    plt.plot (right_arm, E_mm_a_wo_plot, color='olivedrab')
    plt.xlim (-2.5, 2.5)
    plt.ylim (-2.5, 2.5)
    end = time.time()
    plt.savefig ("just-cg-bounds.png", dpi=1200, bbox_inches="tight")

    print (f"Time for finding boundaries is {end-start} seconds")
    x = []
    x.extend (E_mm_n_list)
    x.extend (left_arm)
    x.extend (right_arm)
    y = []
    y.extend (E_mm_a_no_overlap)
    y.extend (E_mm_a_wo_plot)
    y.extend (E_mm_a_wo_plot)

    data = {"E_mm_n": x, "E_mm_a": y}
    df = pd.DataFrame.from_dict (data)
    df.to_csv ("boundaries_total.csv", index=False)

    """
    data = np.vstack ((x, y)).T
    kmeans = KMeans(n_clusters=3)
    kmeans.fit(data)
    labels = kmeans.predict(data)

    # Plot the results
    plt.figure(2)
    plt.scatter(x, y, c=labels, cmap='viridis')
    plt.savefig("division", dpi=1200, bbox_inches="tight")
    """

    """
    y = np.arange (-1.5, 2.5, 0.1)
    x = []

    interval = [-1.9, 2.5]
    E_list   = [0, 0, E_ms_a, E_ms_n]
    for E_mm_a in np.arange (-1.5, 2.5, 0.1):
        E_list [0] = E_mm_a
        interval = [-1.9, 2.5]
        E_mm_n = bisection (E_list, pv, interval, 1000)
        x.append (E_mm_n)
    
    print (f"E_mm_n = {x}")
    # slope, intercept, r_val, p_val, std_err = linregress (x, y)
    # r2 = r_val ** 2
    # print (f"Slope: {slope}")
    # print (f"intercept: {intercept}")
    # print (f"r2: {r2}")

    plt.plot (x, y, marker='o', mec='k')
    plt.xlim (-2.5, 2.5)
    plt.ylim (-2.5, 2.5)
    plt.show ()
    """
