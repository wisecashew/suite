#!~/.conda/envs/data_analysis/bin/python

import numpy as np
import pandas as pd
import copy 
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import scipy.optimize as so
import time
import argparse

# set up the functions that calculate chi_eff

# g is the density of states fraction
g  = 0.5
    
# zij is the partition function 
# fij is the fraction of interactions that are aligned
zmm  = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float64) + (1-g)*np.exp ((-1/T * emmn), dtype = np.float64)
zms  = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float64) + (1-g)*np.exp ((-1/T * emsn), dtype = np.float64)
fmma = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float64) / zmm(emma, emmn, T)
fmsa = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float64) / zms(emsa, emsn, T)

# we now define chi
def chi (emma, emmn, emsa, emsn, pv, T):
    t1 = pv*(fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn) + (1-pv)*emsn
    t2 = pv*(fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn) + (1-pv)*emmn
    return 1/T * (t1 - 0.5 * t2)   

# we define dfij because we need these to write down dchi/dT
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

# now we have defined all of the things that we need to get chi and dchi

# moving on to analysis scripts
# given a chi value, find the c->g/c->c boundaries. 
# this is somewhat of an involved process. 

# as a first step, i will put forward a function that checks if a plot in energy space can show 
# c->g or c->c 

####################################################################################################
####################################################################################################
###
###  SET OF FUNCTION FOR THE C-G BOUNDARIES
###
####################################################################################################
####################################################################################################


def check_temp_interval_for_c_to_g (E_list, pv):
    T   = np.logspace(-2, 2, 100)
    sgn = np.sign( chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 0.01) ) * np.sign(np.max (chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, T) ) )

    if sgn == 1: # if same sign
        return 0

    else: # if opposite sign
        return 1


# Sometimes, what you see is that there might be a region where c->g takes place at a weaker E_mm_n. You need
# to make sure you segregate the weaker E_mm_n region from the region where there is no c->g at weaker E_mm_n

def range_finder_for_c_to_g (E_mm_a_range, E_mm_n_range, E_ms_list, pv):

    no_overlap   = []
    with_overlap = []


    for E_mm_a in E_mm_a_range:
        E_list = [E_mm_a, E_mm_n_range[0 ], E_ms_list[0], E_ms_list[1]]
        if E_mm_a == -1.6740000000000348:
            print (f"E_list = {E_list}")    
        c1 = check_temp_interval_for_c_to_g (E_list, pv)
        E_list = [E_mm_a, E_mm_n_range[-1], E_ms_list[0], E_ms_list[1]]
        if E_mm_a == -1.6740000000000348:
            print (f"E_list = {E_list}")
        c2 = check_temp_interval_for_c_to_g (E_list, pv)

        if E_mm_a == -1.6740000000000348:
            print (f"c1 = {c1}")
            print (f"c2 = {c2}")

        if c1 + c2 == 1:
            no_overlap.append   (E_mm_a)
        else:
            with_overlap.append (E_mm_a)

    return no_overlap, with_overlap


def bisection_for_c_to_g (E_list, pv, interval, itermax):

    for i in range(itermax):

        E_list_l    = copy.copy(E_list) 
        E_list_l[1] = interval[0]
        E_list_r    = copy.copy(E_list)
        E_list_r[1] = interval[1]

        res_l = check_temp_interval_for_c_to_g (E_list_l, pv) 
        res_r = check_temp_interval_for_c_to_g (E_list_r, pv) 

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
            res_m = check_temp_interval_for_c_to_g (E_list_m, pv)

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


def line_search_for_c_to_g (E_list, E_mm_n_range, pv):
    
    line_bools = []
    for E_mm_n in E_mm_n_range:
        E_list[1] = E_mm_n
        line_bools.append ( check_temp_interval_for_c_to_g (E_list, pv) )

    line_bools = np.array (line_bools)
    line_bools = line_bools==0
    if len(line_bools[line_bools]) == 0:
        return [False]
    else:
        E_left  = E_mm_n_range[line_bools][0 ]
        E_right = E_mm_n_range[line_bools][-1]
        return [True, (E_left, E_right)]


def get_slope_top_left_c_to_g (infile1, infile2, outfile):

    df = pd.read_csv (infile1, names=["Emma", "Emmn"], sep=',', engine='python', skiprows=1)

    x = df["Emmn"].values
    y = df["Emma"].values

    X = np.array ([x, y])
    range_of_indices = [0, len(X[0])-1]
    r2_list = []

    for j in range(2):
        L = (range_of_indices[1] - range_of_indices[0]) // 4
        r2_list.clear() 
        for i in range(4):
            m, b, r, p, se = linregress (X[0,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L], X[1,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L])
            r2 = r**2
            r2_list.append (r2)

        elbow_loc = np.argmin (r2_list) 
        range_of_indices = [range_of_indices[0]+elbow_loc * L, range_of_indices[0]+(elbow_loc+1)*L] # np.arange (elbow_loc*L,(elbow_loc+1)*L)

    m1, b1, r1, p1, se1 = linregress (X[0, 0:range_of_indices[0]], X[1, 0:range_of_indices[0]])
    m2, b2, r2, p2, se2 = linregress (X[0, range_of_indices[1]:] , X[1, range_of_indices[1]:] )

    f = open (outfile, 'a')
    f.write (f"Slope of bottom part of the c-g boundary = {m1}, and R-squared = {r1**2}\n")
    f.write (f"Slope of left part of the c-g boundary   = {m2}, and R-squared = {r2**2}\n")


    df = pd.read_csv (infile2, names=["Emma", "Emmn_left", "Emmn_right"], sep=',', engine='python', skiprows=1)

    x = df["Emmn_right"].values
    y = df["Emma"].values 

    m3, b3, r3, p3, se3 = linregress (x, y)    
    f.write (f"Slope of the right part of the c-g boundary = {m3}, and R-squared = {r3**2}\n")
    f.close ()

    return

def generate_boundaries_for_c_to_g (E_ms_a, E_ms_n, pv, outfile, ax):


    E_ms_list = [E_ms_a, E_ms_n]
    E_ms_min  = np.min (E_ms_list)
    E_mm_min  = E_ms_min * 2

    E_mm_a_range = np.arange (E_mm_min+0.01, 2.5, 0.001)
    E_mm_a_range[-1] = 2.5
    E_mm_n_range = np.arange (E_mm_min+0.01, 2.5, 0.001)
    E_mm_n_range[-1] = 2.5

    E_mm_a_no_overlap, E_mm_a_with_overlap = range_finder_for_c_to_g ( E_mm_a_range, E_mm_n_range, E_ms_list, pv)

    E_mm_n_list = []

    # for points with no c->g taking place on the other side, use bisection
    for E_mm_a in E_mm_a_no_overlap:
        E_list = [E_mm_a, 0, E_ms_a, E_ms_n]
        E_mm_n_interval = [E_mm_min+0.01, 2.5]
        E_mm_n_list.append ( bisection_for_c_to_g (E_list, pv, E_mm_n_interval, 100) )

    # we have the top half of the c->g curve
    # let's put it out so I can look at it later 
    data = {"E_mm_a": E_mm_a_no_overlap, "E_mm_n": E_mm_n_list}
    df   = pd.DataFrame.from_dict (data)
    infile1 = "boundaries_top_left_half_c_to_g_EMSA_"+str(E_ms_a)+"_EMSN_"+str(E_ms_n)+".csv"
    df.to_csv (infile1, index=False)

    ax.plot (E_mm_n_list, E_mm_a_no_overlap, c='coral', lw=1, linestyle='-', label="CG/CC")

    # now, let's get the bottom half. 
    left_arm  = [] 
    right_arm = []
    E_mm_a_wo_plot = []

    for E_mm_a in E_mm_a_with_overlap:
        E_list = [E_mm_a, 0, E_ms_a, E_ms_n]
        results = line_search_for_c_to_g (E_list, E_mm_n_range, pv)

        if results[0]:
            left_arm.append  (results[1][0])
            right_arm.append (results[1][1])
            E_mm_a_wo_plot.append (E_mm_a)
        else:
            continue


    ax.plot (left_arm , E_mm_a_wo_plot, c='coral', lw=1, linestyle='-', label="_nolabel_")
    ax.plot (right_arm, E_mm_a_wo_plot, c='coral', lw=1, linestyle='-', label="_nolabel_")

    data = {"E_mm_a": E_mm_a_wo_plot, "E_mm_n_left": left_arm, "E_mm_n_right": right_arm}
    df   = pd.DataFrame.from_dict (data)
    infile2 = "boundaries_bottom_right_half_c_to_g_EMSA_"+str(E_ms_a)+"_EMSN_"+str(E_ms_n)+".csv"
    df.to_csv (infile2, index=False) 

    get_slope_top_left_c_to_g (infile1, infile2, outfile)

    return 


####################################################################################################
####################################################################################################
###
###  END OF FUNCTION FOR THE C-> (G/C) BOUNDARIES
###
####################################################################################################
####################################################################################################


####################################################################################################
####################################################################################################
###
###  SET OF FUNCTION FOR THE G-> (C/G) BOUNDARIES
###
####################################################################################################
####################################################################################################

def check_temp_interval_for_g_to_c (E_list, pv):

    sgn = np.sign( chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 0.01) ) * np.sign(chi (E_list[0], E_list[1], E_list[2], E_list[3], pv, 100))

    if sgn == 1:
        return 0

    else:
        return 1


def bisection_for_g_to_c_horizontal (E_list, pv, interval, itermax):

    for i in range (itermax):

        E_list_l    = copy.copy (E_list)
        E_list_l[1] = interval[0]
        E_list_r    = copy.copy (E_list)
        E_list_r[1] = interval[1]

        res_l = check_temp_interval_for_g_to_c (E_list_l, pv)
        res_r = check_temp_interval_for_g_to_c (E_list_r, pv)

        if res_l + res_r != 1:
            print ("Bad bounds for horizontal sweep!")
            exit  ()
            return 

        else: 
            E_list_m    = copy.copy (E_list)
            E_list_m[1] = np.mean (interval)
            res_m = check_temp_interval_for_g_to_c (E_list_m, pv)

        if res_l + res_m == 1:
            interval[1] = np.mean (interval)

        else:
            interval[0] = np.mean (interval)

        chi_l = lambda T: chi (E_list[0], interval[1], E_list[2], E_list[3], pv, T)
        sol = so.brentq (chi_l, 0.01, 100)        

        if np.abs (100-sol) < 1e-4:
            return interval[1]

        else:
            E_list[1] = interval[1] 


    return interval[1]


def bisection_for_g_to_c_vertical (E_list, pv, interval, itermax):

    for i in range (itermax):

        E_list_l    = copy.copy (E_list)
        E_list_l[0] = interval[0]
        E_list_r    = copy.copy (E_list)
        E_list_r[0] = interval[1]

        res_l = check_temp_interval_for_g_to_c (E_list_l, pv)
        res_r = check_temp_interval_for_g_to_c (E_list_r, pv)

        if res_l + res_r != 1:
            print ("Bad bounds for vertical sweep!")
            exit  ()
            return None

        else:
            E_list_m = copy.copy  (E_list)
            E_list_m[0] = np.mean (interval)
            res_m       = check_temp_interval_for_g_to_c (E_list_m, pv)

        if res_l + res_m == 1:
            interval[1] = np.mean (interval)

        else:
            interval[0] = np.mean (interval)

        chi_l = lambda T: chi(interval[1], E_list[1], E_list[2], E_list[3], pv, T)
        sol   = so.brentq (chi_l, 0.01, 100)

        if np.abs (100-sol) < 1e-4:
            return interval[1]

        else:
            E_list[0] = interval[1]


    return interval[1]


def generate_boundaries_for_g_to_c (E_ms_a, E_ms_n, pv, outfile, ax):

    E_ms_list = [E_ms_a, E_ms_n]
    E_ms_min  = np.min (E_ms_list)
    E_mm_max  = E_ms_min * 2

    # get the E_mm_a's for the horizontal search 
    E_mm_a_range = np.arange (-2.5, E_mm_max-0.01, 0.01)
    E_mm_n_sol   = []

    # get the E_mm_n's for the vertical search 
    E_mm_n_range = np.arange (-2.5, E_mm_max-0.01, 0.01)
    E_mm_a_sol   = []

    E_list = [0, 0, E_ms_a, E_ms_n]

    # perform the vertical search 
    for E_mm_n in E_mm_n_range:
        E_list [1] = E_mm_n
        search_interval = [-3.5, 3.5]
        E_mm_a_to_append = bisection_for_g_to_c_vertical (E_list, pv, search_interval, 1000)
        E_mm_a_sol.append (E_mm_a_to_append)

    m_vert, b_vert, r_vert, p_vert, s_vert = linregress (E_mm_n_range, E_mm_a_sol)

    E_list = [0, 0, E_ms_a, E_ms_n]

    for E_mm_a in E_mm_a_range:
        E_list [0] = E_mm_a
        search_interval = [-3.5, 3.5]
        E_mm_n_to_append = bisection_for_g_to_c_horizontal (E_list, pv, search_interval, 1000)
        E_mm_n_sol.append (E_mm_n_to_append)

    m_hori, b_hori, r_hori, p_hori, s_hori = linregress (E_mm_n_sol, E_mm_a_range)

    ax.plot (E_mm_n_sol, E_mm_a_range, color='steelblue', lw=1, solid_capstyle="round", label="GG/CG")
    ax.plot (E_mm_n_range, E_mm_n_sol, color='steelblue', lw=1, solid_capstyle="round", label="_nolabel_")

    f = open (outfile, 'a')
    f.write (f"Slope of vertical search portion of the g-c boundary   = {m_vert}, and R-squared = {r_vert**2}\n")
    f.write (f"Slope of horizontal search portion of the g-c boundary = {m_hori}, and R-squared = {r_hori**2}\n")
    f.close () 

    return 

####################################################################################################
####################################################################################################
###
###  END OF FUNCTION FOR THE G-> (C/G) BOUNDARIES
###
####################################################################################################
####################################################################################################

parser = argparse.ArgumentParser (description="Create boundary plots and calculate slopes.")
parser.add_argument ("--EMSA", dest='ea', action='store', type=float, help="Value of aligned E_MS energy.")
parser.add_argument ("--EMSN", dest='en', action='store', type=float, help="Value of misaligned E_MS energy.")
args = parser.parse_args ()

if __name__=="__main__":
    
    start = time.time()
    # set up the file which will hold all our results
    outfile = "results_EA_" + str(args.ea) + "_EN_" + str(args.en) + ".out"
    f = open (outfile, 'w')
    f.close ()

    E_ms_n  = args.ea
    E_ms_a  = args.en
    pv      = 1.0

    E_min   = np.min ([E_ms_a, E_ms_n])

    fig = plt.figure (figsize=(4/1.6,3/1.6), constrained_layout=True)
    ax  = plt.axes   ()

    ax.set_ylim (-2.5, 2.5)
    ax.set_xlim (-2.5, 2.5)
    ax.set_xticks ([-2.5, -1, 0, 1, 2.5])
    ax.set_yticks ([-2.5, -1, 0, 1, 2.5])

    ax.plot ( [2*E_min, 2*E_min], [2*E_min, 2.5], c='dimgray', lw=0.5, linestyle='-', label="Primary")
    ax.plot ( [2*E_min, 2.5], [2*E_min, 2*E_min], c='dimgray', lw=0.5, linestyle='-', label="_nolabel_")

    generate_boundaries_for_g_to_c (E_ms_a, E_ms_n, pv, outfile, ax)

    if np.abs(E_ms_n - E_ms_a) > 1e-4:
        generate_boundaries_for_c_to_g (E_ms_a, E_ms_n, pv, outfile, ax)
    else:
        pass

    ax.legend(loc="upper right", prop={'size':5}, frameon=False, fancybox=True)

    fig.savefig ("boundaries_EA_" + str(args.ea) + "_EN_" + str(args.en) + ".png", dpi=1200, bbox_inches="tight")
    end = time.time()

    print (f"Time required for computation is {end-start} seconds.")



