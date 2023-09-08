import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
from scipy.spatial.distance import cdist
import sys
import argparse
import linecache
import mpltern
np.set_printoptions(threshold=sys.maxsize)
import warnings 

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

import argparse 

parser = argparse.ArgumentParser(description="Locate the /c/ritical point on /spin/odal diagram. This program will create two spinodal diagram: one with only shows the edges (ternary) and another one which paints the ternary plot.")
parser.add_argument('--chisc', metavar='chi_sc', dest='chi_sc', type=float, action='store', help='enter A-C exchange parameter.' )
parser.add_argument('--chips', metavar='chi_ps', dest='chi_ps', type=float, action='store', help='enter A-B exchange parameter.' )
parser.add_argument('--chipc', metavar='chi_pc', dest='chi_pc', type=float, action='store', help='enter B-C exchange parameter.' )
parser.add_argument('-N',      metavar='N',      dest='N',      type=int,   action='store', help='degree of polymerization of B.')
parser.add_argument('--ternary', action='store_true', default=False, help='make the output a ternary plot.')
parser.add_argument('--tang_norm', action='store_true', default=False, help='draw normal and tangent at critical point.')
args = parser.parse_args()


def remove_close_rows(array, threshold):

    filtered_array = np.empty ((0,2))
    for i, elem in enumerate(array):
        if i == 0:
            filtered_array = np.vstack((filtered_array, elem))
            continue
        else:
            sieve = (np.linalg.norm(filtered_array - elem, axis=1) < 1e-3).any()
            if sieve:
                continue
            else:
                filtered_array = np.vstack((filtered_array, elem))

    return filtered_array



def crit_condition (N, phi_p, phi_s, chi_sc, chi_ps, chi_pc):

    phi_c = 1-phi_p-phi_s
    t1    = 1/phi_c + 1/(phi_p*N) - 2*chi_pc
    t2    = (1/(phi_c)**2 - 1/phi_s**2)*(1/(phi_c) + 1/(phi_p*N) - 2*chi_pc) + (1/(phi_c) + 1/phi_s - 2*chi_sc)/(phi_c)**2 - 2*(1/phi_c - chi_pc - chi_sc + chi_ps)/(phi_c)**2

    u1    = (1/phi_c + 1/(phi_p*N) - 2*chi_pc)/(phi_c)**2 + (1/(phi_c)**2 - 1/(phi_p**2 * N))*(1/phi_c + 1/phi_s - 2*chi_sc) - 2*(1/phi_c + chi_ps - chi_sc - chi_pc)/phi_c**2
    u2    = 1/phi_c - chi_pc - chi_sc + chi_ps

    return t1*t2 - u1*u2


def find_crit_point (N, chi_sc, chi_ps, chi_pc):

    def send_to_fsolve_r1 (phi_s):
        phi_p_upper = root_up (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (N, phi_p_upper, phi_s, chi_sc, chi_ps, chi_pc)

    def send_to_fsolve_r2 (phi_s):
        phi_p_lower = root_lo (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (N, phi_p_lower, phi_s, chi_sc, chi_ps, chi_pc)

    guesses = np.linspace (0, 1, 10000)
    roots_up   = np.empty ((0,2))
    roots_down = np.empty ((0,2))

    for g in guesses:
        root = fsolve (send_to_fsolve_r1, g)
        
        if abs(send_to_fsolve_r1(root)) < 1e-6:

            if root >= 1 or root <= 0 or np.isnan(root):
                pass
            else:
                r_up  = root_up(root, chi_ps, chi_pc, chi_sc)[0]
                r_tup = np.array([root[0], r_up])
                if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_up):
                    pass

                elif r_tup in roots_up:
                    pass

                else:
                    if len(roots_up) == 0:
                        roots_up = np.vstack ((roots_up,r_tup))
                    else:
                        similarity = (np.linalg.norm(roots_up - r_tup, axis=1) < 1e-3).any ()
                        if similarity:
                            pass
                        else:
                            roots_up = np.vstack ((roots_up,r_tup))   
        else:
            pass

    for g in guesses:
        root = fsolve (send_to_fsolve_r2, g)

        if abs(send_to_fsolve_r2(root)) < 1e-6:

            if root >= 1 or root <= 0 or np.isnan(root):
                pass
            else:
                r_lo = root_lo(root, chi_ps, chi_pc, chi_sc)[0]
                r_tup = np.array([root[0], r_lo])
                if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_lo):
                    pass

                elif r_tup in roots_down:
                    pass

                else: 
                    if len(roots_down) == 0:
                        roots_down = np.vstack ((roots_down,r_tup))
                    else:
                        similarity = (np.linalg.norm(roots_down - r_tup, axis=1) < 1e-3).any ()
                        if similarity:
                            pass
                        else:
                            roots_down = np.vstack ((roots_down,r_tup))

        else:
            pass

    return roots_up, roots_down 




if __name__=="__main__":

    print ("Start painting the spinodal region...")

    chi_sc = args.chi_sc
    chi_ps = args.chi_ps
    chi_pc = args.chi_pc
    N      = args.N


    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4* N * (1 - 2* phi_s * chi_sc + 2 * phi_s ** 2 * chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    (-1 + 2 * phi_s * chi_sc + N * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc:  1 / (2*N * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: 1 - 2 * phi_s * chi_sc + N * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )

    roots_up, roots_down = find_crit_point (N, chi_sc, chi_ps, chi_pc)

    print (f"roots_up   = {roots_up}")
    print (f"roots_down = {roots_down}")


    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    def tangent (ps, pp, vp, cpc, cps, csc):

        print (f"ps = {ps}, pp = {pp}")

        dist_lo  = np.linalg.norm(pp - root_lo (ps, cps, cpc, csc))
        dist_up  = np.linalg.norm(pp - root_up (ps, cps, cpc, csc))

        print (f"lower root distance = {dist_lo}")
        print (f"upper root distance = {dist_up}")

        if dist_lo > dist_up:
            tang_slope = (2 * csc + 2 * cpc * vp - cpc**2 * vp - 2 * cps * vp + 2 * cpc * cps * vp - cps**2 * vp + \
            2 * cpc * csc * vp + 2 * cps * csc * vp - csc**2 * vp + 2 * cpc**2 * ps * vp - \
            4 * cpc * cps * ps * vp + 2 * cps**2 * ps * vp - 4 * cpc * csc * ps * vp - \
            4 * cps * csc * ps * vp + 2 * csc**2 * ps * vp - (-4 * (-1 + 2 * csc * ps - 2 * csc * ps**2) * (-cpc**2 * vp + \
            2 * cpc * cps * vp - cps**2 * vp + 2 * cpc * csc * vp + 2 * cps * csc * vp - csc**2 * vp) - 4 * (2 * csc - 4 * csc * ps) * (-2 * cpc * vp - cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp) + \
            2 * (-2 * csc - 2 * cpc * vp + cpc**2 * vp + 2 * cps * vp - 2 * cpc * cps * vp + cps**2 * vp - 2 * cpc * csc * vp - 2 * cps * csc * vp + csc**2 * vp - 2 * cpc**2 * ps * vp + 4 * cpc * cps * ps * vp - 2 * cps**2 * ps * vp + \
            4 * cpc * csc * ps * vp + 4 * cps * csc * ps * vp - 2 * csc**2 * ps * vp) * (1 - 2 * csc * ps - vp + 2 * cpc * vp - 2 * cpc * ps * vp + cpc**2 * ps * vp + 2 * cps * ps * vp - 2 * cpc * cps * ps * vp + cps**2 * ps * vp - \
            2 * cpc * csc * ps * vp - 2 * cps * csc * ps * vp + csc**2 * ps * vp - cpc**2 * ps**2 * vp + 2 * cpc * cps * ps**2 * vp - cps**2 * ps**2 * vp + 2 * cpc * csc * ps**2 * vp + 2 * cps * csc * ps**2 * vp - csc**2 * ps**2 * vp))/ \
            (2 * np.sqrt(-4 * (-1 + 2 * csc * ps - \
            2 * csc * ps**2) * (-2 * cpc * vp - cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - \
            cps**2 * ps * vp + 2 * cpc * csc * ps * vp + 2 * cps * csc * ps * vp - \
            csc**2 * ps * vp) + (1 - 2 * csc * ps - vp + 2 * cpc * vp - \
            2 * cpc * ps * vp + cpc**2 * ps * vp + 2 * cps * ps * vp - \
            2 * cpc * cps * ps * vp + cps**2 * ps * vp - 2 * cpc * csc * ps * vp - \
            2 * cps * csc * ps * vp + csc**2 * ps * vp - cpc**2 * ps**2 * vp + \
            2 * cpc * cps * ps**2 * vp - cps**2 * ps**2 * vp + 2 * cpc * csc * ps**2 * vp + \
            2 * cps * csc * ps**2 * vp - csc**2 * ps**2 * vp)**2)))/(2 * (-2 * cpc * vp - \
            cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp)) - ((-cpc**2 * vp + 2 * cpc * cps * vp - \
            cps**2 * vp + 2 * cpc * csc * vp + 2 * cps * csc * vp - csc**2 * vp) * (-1 + \
            2 * csc * ps + vp - 2 * cpc * vp + 2 * cpc * ps * vp - cpc**2 * ps * vp - \
            2 * cps * ps * vp + 2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp + cpc**2 * ps**2 * vp - \
            2 * cpc * cps * ps**2 * vp + cps**2 * ps**2 * vp - 2 * cpc * csc * ps**2 * vp - \
            2 * cps * csc * ps**2 * vp + \
            csc**2 * ps**2 * vp - np.sqrt(-4 * (-1 + 2 * csc * ps - \
            2 * csc * ps**2) * (-2 * cpc * vp - cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - \
            cps**2 * ps * vp + 2 * cpc * csc * ps * vp + 2 * cps * csc * ps * vp - \
            csc**2 * ps * vp) + (1 - 2 * csc * ps - vp + 2 * cpc * vp - \
            2 * cpc * ps * vp + cpc**2 * ps * vp + 2 * cps * ps * vp - 2 * cpc * cps * ps * vp +\
            cps**2 * ps * vp - 2 * cpc * csc * ps * vp - 2 * cps * csc * ps * vp + \
            csc**2 * ps * vp - cpc**2 * ps**2 * vp + 2 * cpc * cps * ps**2 * vp - \
            cps**2 * ps**2 * vp + 2 * cpc * csc * ps**2 * vp + 2 * cps * csc * ps**2 * vp - \
            csc**2 * ps**2 * vp)**2)))/(2 * (-2 * cpc * vp - cpc**2 * ps * vp + \
            2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp)**2) 
            print (f"tang_slope = {tang_slope}")

        else:
            tang_slope = (2 * csc + 2 * cpc * vp - cpc**2 * vp - 2 * cps * vp + 2 * cpc * cps * vp - cps**2 * vp + \
            2 * cpc * csc * vp + 2 * cps * csc * vp - csc**2 * vp + 2 * cpc**2 * ps * vp - \
            4 * cpc * cps * ps * vp + 2 * cps**2 * ps * vp - 4 * cpc * csc * ps * vp - \
            4 * cps * csc * ps * vp + \
            2 * csc**2 * ps * vp + (-4 * (-1 + 2 * csc * ps - 2 * csc * ps**2) * (-cpc**2 * vp + \
            2 * cpc * cps * vp - cps**2 * vp + 2 * cpc * csc * vp + 2 * cps * csc * vp - \
            csc**2 * vp) - \
            4 * (2 * csc - 4 * csc * ps) * (-2 * cpc * vp - cpc**2 * ps * vp + \
            2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp) + \
            2 * (-2 * csc - 2 * cpc * vp + cpc**2 * vp + 2 * cps * vp - 2 * cpc * cps * vp + \
            cps**2 * vp - 2 * cpc * csc * vp - 2 * cps * csc * vp + csc**2 * vp - \
            2 * cpc**2 * ps * vp + 4 * cpc * cps * ps * vp - 2 * cps**2 * ps * vp + \
            4 * cpc * csc * ps * vp + 4 * cps * csc * ps * vp - 2 * csc**2 * ps * vp) * (1 - \
            2 * csc * ps - vp + 2 * cpc * vp - 2 * cpc * ps * vp + cpc**2 * ps * vp + \
            2 * cps * ps * vp - 2 * cpc * cps * ps * vp + cps**2 * ps * vp - \
            2 * cpc * csc * ps * vp - 2 * cps * csc * ps * vp + csc**2 * ps * vp - \
            cpc**2 * ps**2 * vp + 2 * cpc * cps * ps**2 * vp - cps**2 * ps**2 * vp + \
            2 * cpc * csc * ps**2 * vp + 2 * cps * csc * ps**2 * vp - \
            csc**2 * ps**2 * vp))/(2 * np.sqrt(-4 * (-1 + 2 * csc * ps - \
            2 * csc * ps**2) * (-2 * cpc * vp - cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - \
            cps**2 * ps * vp + 2 * cpc * csc * ps * vp + 2 * cps * csc * ps * vp - \
            csc**2 * ps * vp) + (1 - 2 * csc * ps - vp + 2 * cpc * vp - \
            2 * cpc * ps * vp + cpc**2 * ps * vp + 2 * cps * ps * vp - \
            2 * cpc * cps * ps * vp + cps**2 * ps * vp - 2 * cpc * csc * ps * vp - \
            2 * cps * csc * ps * vp + csc**2 * ps * vp - cpc**2 * ps**2 * vp + \
            2 * cpc * cps * ps**2 * vp - cps**2 * ps**2 * vp + 2 * cpc * csc * ps**2 * vp + \
            2 * cps * csc * ps**2 * vp - csc**2 * ps**2 * vp)**2)))/(2 * (-2 * cpc * vp - \
            cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp)) - ((-cpc**2 * vp + 2 * cpc * cps * vp - \
            cps**2 * vp + 2 * cpc * csc * vp + 2 * cps * csc * vp - csc**2 * vp) * (-1 + \
            2 * csc * ps + vp - 2 * cpc * vp + 2 * cpc * ps * vp - cpc**2 * ps * vp - \
            2 * cps * ps * vp + 2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp + cpc**2 * ps**2 * vp - \
            2 * cpc * cps * ps**2 * vp + cps**2 * ps**2 * vp - 2 * cpc * csc * ps**2 * vp - \
            2 * cps * csc * ps**2 * vp + \
            csc**2 * ps**2 * vp + np.sqrt (-4 * (-1 + 2 * csc * ps - \
            2 * csc * ps**2) * (-2 * cpc * vp - cpc**2 * ps * vp + 2 * cpc * cps * ps * vp - \
            cps**2 * ps * vp + 2 * cpc * csc * ps * vp + 2 * cps * csc * ps * vp - \
            csc**2 * ps * vp) + (1 - 2 * csc * ps - vp + 2 * cpc * vp - \
            2 * cpc * ps * vp + cpc**2 * ps * vp + 2 * cps * ps * vp - 2 * cpc * cps * ps * vp + \
            cps**2 * ps * vp - 2 * cpc * csc * ps * vp - 2 * cps * csc * ps * vp + \
            csc**2 * ps * vp - cpc**2 * ps**2 * vp + 2 * cpc * cps * ps**2  * vp - \
            cps**2 * ps**2 * vp + 2 * cpc * csc * ps**2 * vp + 2 * cps * csc * ps**2 * vp - \
            csc**2 * ps**2 * vp)**2))) / (2 * (-2 * cpc * vp - cpc**2 * ps * vp + \
            2 * cpc * cps * ps * vp - cps**2 * ps * vp + 2 * cpc * csc * ps * vp + \
            2 * cps * csc * ps * vp - csc**2 * ps * vp)**2)
            print (f"tang_slope = {tang_slope}")

        return tang_slope

    p_s_space = np.arange (0.001, 1-0.001, 0.001)
    p_s = np.repeat (p_s_space, len(p_s_space))

    p_p = np.zeros (p_s.shape)
    for i in range (len(p_s_space)):
        p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

    vals = stab_crit (p_s, p_p, chi_ps, chi_pc, chi_sc)

    to_keep = ~np.isnan(vals)

    vals = vals [to_keep]
    p_s  = p_s  [to_keep]
    p_p  = p_p  [to_keep]

    vmax = np.max (vals)
    vmin = np.min (vals)

    if np.sign (vmax) == np.sign (vmin):
        if np.sign (vmax) >=0:
            vmin = -vmax
        else:
            vmax = -vmin

    if -vmin == vmax:
        print (f"There is no unstable region.")
    

    norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
    cols = cm.bwr (norm (vals))

    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(8,8))
    if args.ternary:
        ax = fig.add_subplot (projection="ternary")
    else:
        ax = plt.axes()

    # Plot the points
    p_c = 1 - p_s - p_p
    crits      = np.vstack ((roots_up, roots_down))
    print (f"crits = \n{crits}")
    threshold  = 1e-3
    crits      = remove_close_rows (crits, threshold)
    print (f"cleaned_crits = \n{crits}")

    if args.ternary:
        ax.scatter (p_s, 1-p_p-p_s, p_p, s=1, color=cols)
        ax.scatter (crits[:,0], 1-crits[:,0]-crits[:,1], crits[:,1], color='forestgreen', edgecolors='forestgreen', s=4)
        ax.scatter (np.mean (crits, axis=0)[0], 1-np.mean(crits, axis=0)[0]-np.mean(crits, axis=0)[1], np.mean(crits,axis=0)[1], color='darkred', edgecolors='coral', s=4)
    else:
        ax.scatter (p_s, p_p, s=1, color=cols)
        ax.scatter (crits[:,0], crits[:,1], color='k', edgecolors='limegreen', s=2)
        ax.scatter (np.mean (crits, axis=0)[0], np.mean(crits, axis=0)[1], color='k', edgecolors='darkred', s=2)
    
    phi_ss = np.linspace (0,1,100)

    
    if args.tang_norm:
        cols = ["coral", "steelblue", "limegreen", "pink", "maroon"]
        for idx,crit in enumerate(crits):    

            slope               = tangent (crit[0], crit[1], N, chi_pc, chi_ps, chi_sc)
            perp_slope          = -1/slope
            tangent_vector      = np.array([1, slope]) / np.sqrt(1+slope**2)
            normal_vector       = np.array([1, perp_slope]) / np.sqrt(1+perp_slope**2)
        
            points_along_tangent = lambda L: np.array([crit[0] + L * tangent_vector[0], crit[1] + L * tangent_vector[1]])
            points_along_normal  = lambda L: np.array([crit[0] + L * normal_vector [0], crit[1] + L * normal_vector[1] ])

            l = np.linspace (-10, 10, 100)
            if args.ternary:
                ax.plot (points_along_tangent(l)[0], 1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_tangent(l)[1], c=cols[idx], lw=0.5)
                ax.plot (points_along_normal(l)[0],  1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_normal(l)[1],  c=cols[idx], lw=0.5)

            else:
                ax.plot (points_along_tangent(l)[0], points_along_tangent(l)[1], c=cols[idx], lw=0.5)
                ax.plot (points_along_normal(l)[0],  points_along_normal(l)[1],  c=cols[idx], lw=0.5)

    if args.ternary:

        ax.set_tlabel('Vol. frac. A')
        ax.set_llabel('Vol. frac. C')
        ax.set_rlabel('Vol. frac. B')
        ax.set_tlim(0, 1)
        ax.set_llim(0, 1)
        ax.set_rlim(0, 1)
        positions = ['tick1', 'tick2']
        for position in positions:
            ax.taxis.set_ticks_position(position)
            ax.laxis.set_ticks_position(position)
            ax.raxis.set_ticks_position(position)

    else:
        ax.set_xlabel ("$\\phi _{S}$")
        ax.set_ylabel ("$\\phi _{P}$")
        ax.set_xlim   (0,1)
        ax.set_ylim   (0,1)

    ax.grid()

    if args.ternary:
        plt.savefig (f"signs_tern-N_{N}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)
    else:
        plt.savefig (f"signs_reg-N_{N}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png",  dpi=1200)

    print ("Completed heat map computation.")

