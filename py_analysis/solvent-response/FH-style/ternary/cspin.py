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

parser = argparse.ArgumentParser(description="Locate the /c/ritical points on the /spin/odal diagram. This program will create one plot, and you can customize what you want on the plot.")
parser.add_argument('--chisc',  metavar='chi_sc',  dest='chi_sc',  type=float,   action='store', help='enter A-C exchange parameter.' )
parser.add_argument('--chips',  metavar='chi_ps',  dest='chi_ps',  type=float,   action='store', help='enter A-B exchange parameter.' )
parser.add_argument('--chipc',  metavar='chi_pc',  dest='chi_pc',  type=float,   action='store', help='enter B-C exchange parameter.' )
parser.add_argument('-vs',      metavar='vs',      dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',      dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',      dest='vp',      type=float,   action='store', help='specific volume of polymer.')
#
parser.add_argument('--dont-calc-crits',     dest='crits',     action='store_false', default=True,  help='Put this in to make sure critical points are not calculated.')
parser.add_argument('--ternary',             dest='ternary',   action='store_true',  default=False, help='make the output a ternary plot.')
parser.add_argument('--draw-edges-spinodal', dest='edges',     action='store_true',  default=False, help='draw the edges of the spinodal.')
parser.add_argument('--tang_norm',           dest='tang_norm', action='store_true',  default=False, help='draw normal and tangent at critical point.')
parser.add_argument('--img-name',            dest='img',       action='store', type=str,  default="None", help='name of the image to be created (default: all of the inputs in the imagename).')
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



def crit_condition (vs, vc, vp, phi_p, phi_s, chi_sc, chi_ps, chi_pc):

    phi_c = 1-phi_p-phi_s
    t1    = 1/(vc*phi_c) + 1/(phi_p*vp) - 2*chi_pc
    t2    = (1/(vc*(phi_c)**2) - 1/(vs*phi_s**2))*(1/(vc*(phi_c)) + 1/(phi_p*vp) - 2*chi_pc) + (1/(vc*phi_c) + 1/(vs*phi_s) - 2*chi_sc)/(vc*(phi_c)**2) - 2*(1/(vc*phi_c) - chi_pc - chi_sc + chi_ps)/(vc*(phi_c)**2)

    u1    = (1/(vc*phi_c) + 1/(phi_p*vp) - 2*chi_pc)/(vc*(phi_c)**2) + (1/(vc*phi_c**2) - 1/(phi_p**2 * vp))*(1/(vc*phi_c) + 1/(vs*phi_s) - 2*chi_sc) - 2*(1/(vc*phi_c) + chi_ps - chi_sc - chi_pc)/(vc*phi_c**2)
    u2    = 1/(vc*phi_c) - chi_pc - chi_sc + chi_ps

    return t1*t2 - u1*u2


def find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc):

    def send_to_fsolve_r1 (phi_s):
        phi_p_upper = root_up (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (vs, vc, vp, phi_p_upper, phi_s, chi_sc, chi_ps, chi_pc)

    def send_to_fsolve_r2 (phi_s):
        phi_p_lower = root_lo (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (vs, vc, vp, phi_p_lower, phi_s, chi_sc, chi_ps, chi_pc)

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

    print ("Revving up the program...", flush=True)

    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(8,8))
    if args.ternary:
        ax = fig.add_subplot (projection="ternary")
    else:
        ax = plt.axes()

    chi_sc = args.chi_sc
    chi_ps = args.chi_ps
    chi_pc = args.chi_pc
    vs     = args.vs
    vc     = args.vc
    vp     = args.vp


    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc)) + (vp - 2*phi_s*vp *vs *chi_ps + vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc) ) ) )**2 # -4* vb * (- vc + phi_s * vc - phi_s * vs + 2*phi_s*vc*vs*chi_sc - 2*phi_s**2*vc*vs*chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    # (-1 + 2 * phi_s * chi_sc + vb * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))  # 1 / (2*vb * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp * (2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc) ) ) # 1 - 2 * phi_s * chi_sc + vb * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )

    def stab_crit (p_s, p_p, c_ps, c_pc, c_sc):
        return (1/(vp*p_p) + 1/(vc*(1-p_s - p_p)) - 2 * c_pc) * (1/(vs*p_s) + 1/(vc*(1-p_s - p_p)) - 2 * c_sc) - (1/(vc*(1-p_s-p_p)) + c_ps - c_pc - c_sc) ** 2

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


    def tangent2 (vs, vc, vp, ps, pp, cpc, cps, csc):

        # print (f"ps = {ps}, pp = {pp}")
        # print (f"vs = {vs}, vc = {vc}, vp = {vp}, cpc = {cpc}, cps = {cps}, csc = {csc}.")

        dist_lo  = np.linalg.norm(pp - root_lo (ps, cps, cpc, csc))
        dist_up  = np.linalg.norm(pp - root_up (ps, cps, cpc, csc))

        # print (f"lower root distance = {dist_lo}")
        # print (f"upper root distance = {dist_up}")

        if dist_lo > dist_up:
            tang_slope = (-((2 * cpc + ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (2 * vc * vp * cpc - \
            vc * vp * vs * cpc**2 + 2 * ps * vc * vp * vs * cpc**2 - \
            2 * vp * vs * cps + 2 * vc * vp * vs * cpc * cps - \
            4 * ps * vc * vp * vs * cpc * cps - vc * vp * vs * cps**2 + \
            2 * ps * vc * vp * vs * cps**2 + 2 * vc * vs * csc + \
            2 * vc * vp * vs * cpc * csc - \
            4 * ps * vc * vp * vs * cpc * csc + \
            2 * vc * vp * vs * cps * csc - \
            4 * ps * vc * vp * vs * cps * csc - vc * vp * vs * csc**2 + \
            2 * ps * vc * vp * vs * csc**2 + (-4 * vc * vp * vs * (cpc**2 + \
            (cps - csc)**2 - \
            2 * cpc * (cps + csc)) * (ps * vs + (-1 + \
            ps) * vc * (-1 + 2 * ps * vs * csc)) - \
            4 * vc * vp * (2 * cpc + ps * vs * cpc**2 + \
            ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (vs + \
            vc * (-1 + 2 * (-1 + 2 * ps) * vs * csc)) + \
            2 * (vc + vp * (-1 + 2 * ps * vs * cps) - \
            2 * ps * vc * vs * csc - (-1 + ps) * vc * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))) * (2 * vp * vs * \
            cps - 2 * vc * vs * csc + \
            vc * vp * ((vs - 2 * ps * vs) * cpc**2 - (-1 + \
            2 * ps) * vs * (cps - csc)**2 + cpc * (-2 + \
            2 * (-1 + 2 * ps) * vs * (cps + csc)))))/(2 * \
            np.sqrt(-4 * vc * vp * (2 * cpc + ps * vs * cpc**2 + \
            ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (ps * vs + (-1 \
            + ps) * vc * (-1 + 2 * ps * vs * csc)) + (vp - 2 * ps * vp * vs * cps + \
            vc * (-1 + \
            2 * ps * vs * csc + (-1 + ps) * vp * (2 * cpc + \
            ps * vs * cpc**2 + \
            ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))))**2)))) + \
            vs * (cpc**2 + (cps - csc)**2 - \
            2 * cpc * (cps + csc)) * (-vc + vp - \
            2 * vc * vp * cpc + 2 * ps * vc * vp * cpc - \
            ps * vc * vp * vs * cpc**2 + ps**2 * vc * vp * vs * cpc**2 - \
            2 * ps * vp * vs * cps + 2 * ps * vc * vp * vs * cpc * cps - \
            2 * ps**2 * vc * vp * vs * cpc * cps - ps * vc * vp * vs * cps**2 + \
            ps**2 * vc * vp * vs * cps**2 + 2 * ps * vc * vs * csc + \
            2 * ps * vc * vp * vs * cpc * csc - \
            2 * ps**2 * vc * vp * vs * cpc * csc + \
            2 * ps * vc * vp * vs * cps * csc - \
            2 * ps**2 * vc * vp * vs * cps * csc - ps * vc * vp * vs * csc**2 + \
            ps**2 * vc * vp * vs * csc**2 + np.sqrt(-4 * vc * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (ps * vs + (-1 + \
            ps) * vc * (-1 + 2 * ps * vs * csc)) + (vp - \
            2 * ps * vp * vs * cps + \
            vc * (-1 + \
            2 * ps * vs * csc + (-1 + ps) * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))))**2)))/(2 * vc \
            * vp * (2 * cpc + ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))**2)
            print (f"tang_slope = {tang_slope}")

        else:
            tang_slope = (-((2 * cpc + ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (2 * vc * vp * cpc - \
            vc * vp * vs * cpc**2 + 2 * ps * vc * vp * vs * cpc**2 - \
            2 * vp * vs * cps + 2 * vc * vp * vs * cpc * cps - \
            4 * ps * vc * vp * vs * cpc * cps - vc * vp * vs * cps**2 + \
            2 * ps * vc * vp * vs * cps**2 + 2 * vc * vs * csc + \
            2 * vc * vp * vs * cpc * csc - \
            4 * ps * vc * vp * vs * cpc * csc + \
            2 * vc * vp * vs * cps * csc - \
            4 * ps * vc * vp * vs * cps * csc - vc * vp * vs * csc**2 + \
            2 * ps * vc * vp * vs * csc**2 + (4 * vc * vp * vs * (cpc**2 + \
            (cps - csc)**2 - \
            2 * cpc * (cps + csc)) * (ps * vs + (-1 + \
            ps) * vc * (-1 + 2 * ps * vs * csc)) + \
            4 * vc * vp * (2 * cpc + ps * vs * cpc**2 + 
            ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (vs + \
            vc * (-1 + 2 * (-1 + 2 * ps) * vs * csc)) - \
            2 * (vc + vp * (-1 + 2 * ps * vs * cps) - 
            2 * ps * vc * vs * csc - (-1 + ps) * vc * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - 
            2 * ps * vs * cpc * (cps + csc))) * (2 * vp * vs * \
            cps - 2 * vc * vs * csc + \
            vc * vp * ((vs - 2 * ps * vs) * cpc**2 - (-1 + \
            2 * ps) * vs * (cps - csc) **2 + cpc * (-2 + 
            2 * (-1 + 2 * ps) * vs * (cps + csc))))) / (2 * \
            np.sqrt(-4 * vc * vp * (2 * cpc + ps * vs * cpc**2 + \
            ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (ps * vs + (-1 \
            + ps) * vc * (-1 + 2 * ps * vs * csc)) + (vp - 2 * ps * vp * vs * cps + \
            vc * (-1 + 2 * ps * vs * csc + (-1 + ps) * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))))**2)))) + \
            vs * (cpc**2 + (cps - csc)**2 - \
            2 * cpc * (cps + csc)) * (-vc + vp - \
            2 * vc * vp * cpc + 2 * ps * vc * vp * cpc - \
            ps * vc * vp * vs * cpc**2 + ps**2 * vc * vp * vs * cpc**2 - \
            2 * ps * vp * vs * cps + 2 * ps * vc * vp * vs * cpc * cps - \
            2 * ps**2 * vc * vp * vs * cpc * cps - ps * vc * vp * vs * cps**2 + \
            ps**2 * vc * vp * vs * cps**2 + 2 * ps * vc * vs * csc + \
            2 * ps * vc * vp * vs * cpc * csc - \
            2 * ps**2 * vc * vp * vs * cpc * csc + \
            2 * ps * vc * vp * vs * cps * csc - \
            2 * ps**2 * vc * vp * vs * cps * csc - ps * vc * vp * vs * csc**2 + \
            ps**2 * vc * vp * vs * csc**2 - np.sqrt(-4 * vc * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc)) * (ps * vs + (-1 + \
            ps) * vc * (-1 + 2 * ps * vs * csc)) + (vp - \
            2 * ps * vp * vs * cps + \
            vc * (-1 + \
            2 * ps * vs * csc + (-1 + ps) * vp * (2 * cpc + \
            ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))))**2)))/(2 * vc * \
            vp * (2 * cpc + ps * vs * cpc**2 + ps * vs * (cps - csc)**2 - \
            2 * ps * vs * cpc * (cps + csc))**2)
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

    if len(vals) == 0:
        print (f"There will be no critical points and no spinodal region.")

    vmax = np.max (vals)
    vmin = np.min (vals)
    norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
    cols = cm.bwr (norm (vals))

    if np.sign (vmax) == np.sign (vmin):
        if np.sign (vmax) >=0:
            vmin = -vmax
        else:
            vmax = -vmin

    if -vmin == vmax:
        print (f"There is no unstable region.")
    

    if args.crits:
        roots_up, roots_down = find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc)
        crits      = np.vstack ((roots_up, roots_down))
        threshold  = 1e-6
        crits      = remove_close_rows (crits, threshold)
        print (f"cleaned_crits = \n{crits}")
        
    else:
        print (f"We won't be calculating critical points.")

    # Plot the points
    p_c = 1 - p_s - p_p

    if args.ternary:
        ax.scatter (p_s, 1-p_p-p_s, p_p, s=1, color=cols)
        if args.crits:
            ax.scatter (crits[:,0], 1-crits[:,0]-crits[:,1], crits[:,1], color='darkred', edgecolors='darkred', s=4, zorder=15)
            ax.scatter (np.mean (crits, axis=0)[0], 1-np.mean(crits, axis=0)[0]-np.mean(crits, axis=0)[1], np.mean(crits,axis=0)[1], color='darkred', edgecolors='darkred', s=4, zorder=15)
        else: 
            pass
       
        if args.edges:
            meshsize            = 1000
            phi_s               = np.linspace (0.001, 1-0.001, meshsize*10)
            chi_ps              = args.chi_ps
            chi_pc              = args.chi_pc
            chi_sc              = args.chi_sc

            r1 = root_up (phi_s, chi_ps, chi_pc, chi_sc)
            r2 = root_lo (phi_s, chi_ps, chi_pc, chi_sc)

            to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
            r1 = r1[to_keep_1]

            to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
            r2 = r2[to_keep_2]

            # Plot the points
            ax.scatter(phi_s[to_keep_1], 1-phi_s[to_keep_1]-r1, r1, color='springgreen', s=1)
            ax.scatter(phi_s[to_keep_2], 1-phi_s[to_keep_2]-r2, r2, color='darkturquoise', s=1)
            

    else:
        ax.scatter (p_s, p_p, s=1, color=cols)
        if args.crits:
            ax.scatter (crits[:,0], crits[:,1], color='darkred', edgecolors='darkred', s=4, zorder=15)
            ax.scatter (np.mean (crits, axis=0)[0], np.mean(crits, axis=0)[1], color='darkred', edgecolors='darkred', s=4, zorder=15)
        else:
            pass


        if args.edges:
            meshsize            = 1000
            phi_s               = np.linspace (0.001, 1-0.001, meshsize*100)
            chi_ps              = args.chi_ps
            chi_pc              = args.chi_pc
            chi_sc              = args.chi_sc

            r1 = root_up (phi_s, chi_ps, chi_pc, chi_sc)
            r2 = root_lo (phi_s, chi_ps, chi_pc, chi_sc)

            to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
            r1 = r1[to_keep_1]

            to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
            r2 = r2[to_keep_2]

            # Plot the points
            ax.scatter(phi_s[to_keep_1], r1, color='springgreen', s=1)
            ax.scatter(phi_s[to_keep_2], r2, color='darkturquoise', s=1)


    if args.tang_norm and args.crits:
        cols = ["steelblue", "limegreen", "pink", "maroon"]
        for idx,crit in enumerate(crits):    
            idx = idx%len(crits)
            slope               = tangent2 (vs, vc, vp, crit[0], crit[1], chi_pc, chi_ps, chi_sc)
            perp_slope          = -1/slope
            tangent_vector      = np.array([1, slope]) / np.sqrt(1+slope**2)
            normal_vector       = np.array([1, perp_slope]) / np.sqrt(1+perp_slope**2)
        
            points_along_tangent = lambda L: np.array([crit[0] + L * tangent_vector[0], crit[1] + L * tangent_vector[1]])
            points_along_normal  = lambda L: np.array([crit[0] + L * normal_vector [0], crit[1] + L * normal_vector[1] ])

            l = np.linspace (-10, 10, 100)
            if args.ternary:
                ax.plot (points_along_tangent(l)[0], 1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_tangent(l)[1], c=cols[idx], lw=1)
                ax.plot (points_along_normal(l)[0],  1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_normal(l)[1],  c=cols[idx], lw=1)

            else:
                ax.plot (points_along_tangent(l)[0], points_along_tangent(l)[1], c=cols[idx], lw=1)
                ax.plot (points_along_normal(l)[0],  points_along_normal(l)[1],  c=cols[idx], lw=1)

    if args.ternary:

        ax.set_tlabel ("$\\phi _{S}$") # ('Vol. frac. A')
        ax.set_llabel ("$\\phi _{C}$") # ('Vol. frac. C')
        ax.set_rlabel ("$\\phi _{P}$") # ('Vol. frac. B')
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

    if args.img != "None":
        if "." in args.img:
            plt.savefig (args.img+".png", dpi=1200)
        else:
            plt.savefig (args.img, dpi=1200)
    elif args.ternary:
        plt.savefig (f"signs_tern-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)
    else:
        plt.savefig (f"signs_reg-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.png", dpi=1200)

    print ("Completed heat map computation.")

