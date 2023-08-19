import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
from scipy.optimize import brentq
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import sys
import argparse
import linecache
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
args = parser.parse_args()


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
                if r_up >= 1 or r_up <= 0 or np.isnan(r_up):
                    pass

                elif r_tup in roots_up:
                    pass

                else:
                    if len(roots_up) == 0:
                        roots_up = np.vstack ((roots_up,r_tup))
                    else:
                        # print(roots_up - r_tup)
                        similarity = (np.linalg.norm(roots_up - r_tup, axis=1) < 1e-3).any ()
                        if similarity:
                            pass
                        else:
                            roots_up = np.vstack ((roots_up,r_tup))   
        else:
            pass

    for g in guesses:
        root = fsolve (send_to_fsolve_r2, g)
        # print (root)

        if abs(send_to_fsolve_r2(root)) < 1e-6:

            if root >= 1 or root <= 0 or np.isnan(root):
                pass
            else:
                r_lo = root_lo(root, chi_ps, chi_pc, chi_sc)[0]
                r_tup = np.array([root[0], r_lo])
                if r_lo > 1 or r_lo < 0 or np.isnan(r_lo):
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


    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    roots_up, roots_down = find_crit_point (N, chi_sc, chi_ps, chi_pc)
    crit_points = np.vstack((roots_up, roots_down))
    print (f"crits = \n{crit_points}")
    mesh = 200
    phi_b = np.linspace (0.001, 0.999, mesh)

    phi_b = np.repeat (phi_b, mesh)
    phi_a = np.zeros  (phi_b.shape)

    for i in range (mesh):
        upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
        phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

    # only keep stuff which is outside the spinodal
    to_keep = stab_crit (phi_a, phi_b, chi_ps, chi_pc, chi_sc) > 0

    phi_b   = phi_b [to_keep]
    phi_a   = phi_a [to_keep]

    # now start splitting up phi_a, phi_b
    print ("Cutting up the space...")
    phis    = np.vstack((phi_a, phi_b)).T

    center         = np.mean (crit_points, axis=0)[:2]
    central_axis   = (crit_points[0,:2]-center)/np.linalg.norm (crit_points[0,:2]-center)

    # find those ABOVE axis
    direction      = (phis - center) / np.linalg.norm(phis-center, axis=1)[:, np.newaxis]
    clock          = np.cross (central_axis, direction)
    phi_upper      = phis[clock > 0]
    phi_lower      = phis[clock < 0]

    plt.scatter (phi_upper[:,0], phi_upper[:,1], c='salmon', s=0.1)
    plt.scatter (phi_lower[:,0], phi_lower[:,1], c='steelblue', s=0.1)
    plt.plot (crit_points[:,0], crit_points[:,1], c='red', lw=1)
    plt.savefig ("phis.png", dpi=1200, bbox_inches="tight")
    print ("Plotted!")
   
