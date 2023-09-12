import numpy as np
import numba as nb
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
from scipy.optimize import root
import scipy.spatial.distance
from scipy.spatial.distance import cdist
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import sys
import os
import argparse
import time
import warnings
import linecache
import itertools

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--no-rtw', dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--mesh',   metavar='mesh',   dest='mesh', type=int, action='store', help='enter mesh fineness.')
parser.add_argument('--skelfile', dest='skelfile', type=str, action='store', help="name of file where we are going to dump the skeleton of the binodal (the default contains all the information provided above).", default="None")
args = parser.parse_args()

######
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    if args.nrtw:
        return f"beep.\n"
    else:
        return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
######

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

######

######

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


######


#########################################
def perform_sweep (phi_b, mesh, chi_ps, chi_pc, chi_sc, crit_point, phi_a_edge, phi_b_edge, center):

    print (f"pid = {os.getpid()}.", flush=True)
    phi_b = np.repeat (phi_b, mesh)
    phi_a = np.zeros  (phi_b.shape)

    for i in range (mesh):
        upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
        # phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)
        phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

    # only keep stuff which is outside the spinodal
    to_keep = stab_crit (phi_a, phi_b, chi_ps, chi_pc, chi_sc) >= 0

    phi_b   = phi_b [to_keep]
    phi_a   = phi_a [to_keep]

    phi_b   = np.hstack((phi_b, phi_b_edge))
    phi_a   = np.hstack((phi_a, phi_a_edge))

    # now start splitting up phi_a, phi_b
    phis    = np.vstack((phi_a, phi_b)).T

    central_axis   = (crit_point-center) / np.linalg.norm (crit_point-center)

    # find those ABOVE axis
    direction      = (phis - center) / np.linalg.norm(phis-center, axis=1)[:, np.newaxis]
    clock          = np.cross (central_axis, direction)
    phi_upper      = phis[clock > 0]
    phi_lower      = phis[clock < 0]

    # now split phi_a, phi_b such that they are on either side of the critical line

    chem_pot_a_upper = mu_a (phi_upper[:,0], phi_upper[:,1])
    chem_pot_b_upper = mu_b (phi_upper[:,0], phi_upper[:,1])
    chem_pot_c_upper = mu_c (phi_upper[:,0], phi_upper[:,1])

    masks = np.isinf(chem_pot_a_upper) | np.isnan(chem_pot_a_upper) | np.isinf(chem_pot_b_upper) | np.isnan(chem_pot_b_upper) | np.isinf (chem_pot_c_upper) | np.isnan (chem_pot_c_upper)

    # print ("mu a upper is inf = ",np.isinf(chem_pot_a_upper).any())
    chem_pot_a_upper = chem_pot_a_upper [~masks]
    # print ("mu a upper is inf = ",np.isinf(chem_pot_a_upper).any())
    chem_pot_b_upper = chem_pot_b_upper [~masks]
    chem_pot_c_upper = chem_pot_c_upper [~masks]
    phi_upper = phi_upper[~masks]

    chem_pot_a_lower = mu_a (phi_lower[:,0], phi_lower[:,1])
    chem_pot_b_lower = mu_b (phi_lower[:,0], phi_lower[:,1])
    chem_pot_c_lower = mu_c (phi_lower[:,0], phi_lower[:,1])

    masks = np.isinf(chem_pot_a_lower) | np.isnan(chem_pot_a_lower) | np.isinf(chem_pot_b_lower) | np.isnan(chem_pot_b_lower) | np.isinf (chem_pot_c_lower) | np.isnan (chem_pot_c_lower)

    chem_pot_a_lower = chem_pot_a_lower [~masks]
    chem_pot_b_lower = chem_pot_b_lower [~masks]
    chem_pot_c_lower = chem_pot_c_lower [~masks]
    phi_lower = phi_lower[~masks]

    mu_upper   = np.array ([chem_pot_a_upper, chem_pot_b_upper, chem_pot_c_upper]).T
    mu_lower   = np.array ([chem_pot_a_lower, chem_pot_b_lower, chem_pot_c_lower]).T

    distances       = np.linalg.norm (mu_upper[:, np.newaxis] - mu_lower, axis=2)
    print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
    closest_indices = np.argmin(distances, axis=1)

    phi_lower       = phi_lower[closest_indices]
    mu_lower        = mu_lower [closest_indices]
    min_distances   = distances[np.arange(len(mu_upper)), closest_indices]

    return [phi_upper, phi_lower, min_distances]

#########################################

if __name__=="__main__":

    start = time.time()

    chi_sc = args.chi_sc
    chi_ps = args.chi_ps
    chi_pc = args.chi_pc

    vs = args.vs; # va = args.vs
    vp = args.vp; # vb = args.vp
    vc = args.vc; # vc = args.vc


    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc)) + (vp - 2*phi_s*vp *vs *chi_ps + vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc) ) ) )**2

    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp * (2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )

    #######################################
    # find the edges of the binodal curve
    meshsize            = 1000
    phi_a               = np.linspace (0.001, 1-0.001, meshsize*10)
    chi_ps              = args.chi_ps
    chi_pc              = args.chi_pc
    chi_sc              = args.chi_sc

    r1 = root_up (phi_a, chi_ps, chi_pc, chi_sc)
    r2 = root_lo (phi_a, chi_ps, chi_pc, chi_sc)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
    r2 = r2[to_keep_2]

    phi_a_edge_1 = phi_a[to_keep_1]
    phi_a_edge_2 = phi_a[to_keep_2]

    phi_a_edge = np.hstack((phi_a_edge_1,phi_a_edge_2))
    phi_b_edge = np.hstack((r1,r2))
    if len(phi_a_edge) == 0:
        print ("There is no critical region. No binodals will be found. Exiting...")
        exit ()

    ########################################

    roots_up, roots_down = find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc)
    crits = np.vstack((roots_up, roots_down))
    crits = remove_close_rows (crits, 1e-6)

    # if len(crits) > 2:
    #     print ("Number of crit points is greater than 2. The binodal calculation is going to be wacky and unstable. Exiting computation...", flush=True)
    #     exit  ()

    def stab_crit (p_s, p_p, c_ps, c_pc, c_sc):
        return (1/(vp*p_p) + 1/(vc*(1-p_s - p_p)) - 2 * c_pc) * (1/(vs*p_s) + 1/(vc*(1-p_s - p_p)) - 2 * c_sc) - (1/(vc*(1-p_s-p_p)) + c_ps - c_pc - c_sc) ** 2

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
            print (f"tang_slope = {tang_slope}", flush=True)

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
            print (f"tang_slope = {tang_slope}", flush=True)

        return tang_slope

    ##########################

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - vs/vp * phi_b - vs/vc * (1-phi_a-phi_b) + vs * (phi_b**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_sc + phi_b * (1-phi_a-phi_b) * (chi_ps + chi_sc - chi_pc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vp/vs * phi_a - vp/vc * (1-phi_a-phi_b) + vp * (phi_a**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_pc + phi_a * (1-phi_a-phi_b) * (chi_ps + chi_pc - chi_sc) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/vs * phi_a - vc/vp * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )

    mesh  = args.mesh

    if args.skelfile == "None":
        skelfile = f"bin-vs_{vs}-vc_{vc}-vp_{vp}-chisc_{chi_sc}-chips_{chi_ps}-chipc_{chi_pc}.skelfile"
    else:
        skelfile = args.skelfile
    f = open (skelfile, 'w')
    f.write ("{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}\n"\
    .format ("dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2")) 

    phi_b_spin_min = 0.01
    phi_b_spin_max = 0.8

    phi_b_list = [np.linspace(phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1, mesh), \
          np.linspace(phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3, mesh), \
          np.linspace(phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5, mesh), \
          np.linspace(phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9, mesh), \
          np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh), \
          np.logspace(-20, np.log10(0.9), mesh)]

    print (f"crits = {crits}", flush=True)
    for crit in crits:
        print (f"@ crit point = {crit}...", flush=True)
        tangent_to_crit = tangent2  (vs, vc, vp, crit[0], crit[1], chi_ps, chi_pc, chi_sc)
        normal_slope    = -1/tangent_to_crit
        if len(crits) == 1:
            center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit
        elif len(crits) == 2:
            center  = np.mean(crits, axis=0)
        else:
            center  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2) + crit

        f.write  ("{:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f}\n"\
        .format(0, crit[0], crit[1], 1-crit[0]-crit[1], crit[0], crit[1], 1-crit[0]-crit[1]))


        print (f"# of iterations = {len(phi_b_list)}...", flush=True)
        for idx,b_list in enumerate(phi_b_list):
            print (f"\tAt iteration {idx}...", flush=True)
            results = perform_sweep (b_list, mesh, chi_ps, chi_pc, chi_sc, crit, phi_a_edge, phi_b_edge, center)

            if len(results[0]) == 0:
                print (f"\tNo critical condition in process {idx}.", flush=True)
                continue
            for i in range(len (results[0]) ):
                f.write  ("{:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f}\n"\
        .format(results[2][i], results[0][i,0], results[0][i,1], 1-results[0][i,0]-results[0][i,1], results[1][i,0], results[1][i,1], 1-results[1][i,0]-results[1][i,1]))

    f.close ()

    stop = time.time ()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)


