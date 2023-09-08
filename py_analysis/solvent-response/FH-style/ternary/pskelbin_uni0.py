import numpy as np
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
import argparse
import time
import warnings
import linecache
import multiprocessing as mp
import itertools

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a  memory-intensive computation.")
parser.add_argument('--chisc', metavar='chi_sc', dest='chi_sc', type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips', metavar='chi_ps', dest='chi_ps', type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc', metavar='chi_pc', dest='chi_pc', type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('--mesh', metavar='mesh', dest='mesh', type=int, action='store', help='enter mesh fineness.')
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
parser.add_argument('--skelfile', dest='skelfile', type=str, action='store', help="name of file where we are going to dump the skeleton of the binodal.")
args = parser.parse_args()

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

def process_sweep(args):
    try:
        return perform_sweep(*args)
    except Exception as e:
        traceback.print_exc()
        return []

######

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



######

def get_limits (N):

    discriminant = lambda phi_a, chi_ab, chi_bc, chi_ac: -4* N * (1 - 2* phi_a * chi_ac + 2 * phi_a ** 2 * chi_ac) * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) **2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) + \
    (-1 + 2 * phi_a * chi_ac + N * (1 - 2*chi_bc - phi_a * (chi_ab ** 2 + chi_ac **2 - 2*chi_ac*chi_bc + (chi_bc -2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) + phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) ) ** 2

    denom  = lambda phi_a, chi_ab, chi_bc, chi_ac:  1 / (2*N * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2*chi_ab * (chi_ac + chi_bc) ) ) )

    prefac = lambda phi_a, chi_ab, chi_bc, chi_ac: 1 - 2 * phi_a * chi_ac + N * ( -1 + 2 * chi_bc + phi_a * (chi_ab ** 2 + chi_ac ** 2 - 2*chi_ac * chi_bc + (chi_bc - 2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) - phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) )

    root_up  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) + np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )
    root_lo  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) - np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )

    meshsize            = 1000
    phi_a               = np.logspace (-6, np.log10(1-0.001), meshsize*10)
    chi_ab              = args.chi_ab
    chi_bc              = args.chi_bc
    chi_ac              = args.chi_ac

    r1 = root_up (phi_a, chi_ab, chi_bc, chi_ac)
    r2 = root_lo (phi_a, chi_ab, chi_bc, chi_ac)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0) * (r1+phi_a<=1)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0) * (r2+phi_a<=1)
    r2 = r2[to_keep_2]

    if len(r1) == 0 or len(r2) == 0:
        return None, None

    else:
        curve_x = np.hstack((phi_a[to_keep_1], phi_a[to_keep_2]))
        curve_y = np.hstack((r1, r2))
        edges_x = [np.min(curve_x), np.max (curve_x)]
        edges_y = [np.min(curve_y), np.max (curve_y)]
        return edges_x, edges_y

#########################################

def perform_sweep (phi_b, mesh, chi_ab, chi_bc, chi_ac, crit_point, normal_to_crit, phi_a_edge, phi_b_edge):

    phi_b = np.repeat (phi_b, mesh)
    phi_a = np.zeros  (phi_b.shape)

    for i in range (mesh):
        upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
        # phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)
        phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

    # only keep stuff which is outside the spinodal
    to_keep = stab_crit (phi_a, phi_b, chi_ab, chi_bc, chi_ac) >= -0.1

    phi_b   = phi_b [to_keep]
    phi_a   = phi_a [to_keep]

    phi_b   = np.hstack((phi_b, phi_b_edge))
    phi_a   = np.hstack((phi_a, phi_a_edge))

    # now start splitting up phi_a, phi_b
    phis    = np.vstack((phi_a, phi_b)).T
    # print (phis.shape)

    center         = crit_point
    central_axis   = (crit_point + normal_to_crit) / np.linalg.norm (crit_point + normal_to_crit)

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
    try:
        closest_indices = np.argmin(distances, axis=1)
        # print (phi_b)
        # exit ()
    except ValueError:
        print ("This is the spot.")
        print (phi_b)
        exit ()

    phi_lower       = phi_lower[closest_indices]
    mu_lower        = mu_lower [closest_indices]
    min_distances   = distances[np.arange(len(mu_upper)), closest_indices]

    return [phi_upper, phi_lower, min_distances]

#########################################

if __name__=="__main__":

    start = time.time()

    N      = args.N
    chi_ac = args.chi_ac
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc

    va = 1
    vb = N
    vc = 1

    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4* N * (1 - 2* phi_s * chi_sc + 2 * phi_s ** 2 * chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    (-1 + 2 * phi_s * chi_sc + N * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom    = lambda phi_s, chi_ps, chi_pc, chi_sc:  1 / (2*N * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac   = lambda phi_s, chi_ps, chi_pc, chi_sc: 1 - 2 * phi_s * chi_sc + N * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )


    #######################################
    # find the edges of the binodal curve
    meshsize            = 1000
    phi_a               = np.linspace (0.001, 1-0.001, meshsize*10)
    chi_ab              = args.chi_ab
    chi_bc              = args.chi_bc
    chi_ac              = args.chi_ac

    r1 = root_up (phi_a, chi_ab, chi_bc, chi_ac)
    r2 = root_lo (phi_a, chi_ab, chi_bc, chi_ac)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
    r2 = r2[to_keep_2]

    phi_a_edge_1 = phi_a[to_keep_1]
    phi_a_edge_2 = phi_a[to_keep_2]

    phi_a_edge = np.hstack((phi_a_edge_1,phi_a_edge_2))
    phi_b_edge = np.hstack((r1,r2))

    ########################################

    roots_up, roots_down = find_crit_point (N, chi_ac, chi_ab, chi_bc)
    crits = np.vstack((roots_up, roots_down))
    crits = remove_close_rows (crits, 1e-3)

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    def tangent (ps, pp, vp, cpc, cps, csc):

        dist_lo  = np.linalg.norm(pp - root_lo (ps, cps, cpc, csc))
        dist_up  = np.linalg.norm(pp - root_up (ps, cps, cpc, csc))

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

        return tang_slope

    ##########################

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    mesh  = args.mesh

    edges_x, edges_y = get_limits (N)
    f = open (args.skelfile, 'w')
    f.write ("{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}\n"\
    .format ("dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2")) 

    for crit in crits:
        print (f"@ crit point = {crit}...", flush=True)
        tangent_to_crit = tangent  (crit[0], crit[1], vb, chi_ab, chi_bc, chi_ac)
        normal_slope    = -1/tangent_to_crit
        normal_to_crit  = np.array ([1, normal_slope]) / np.sqrt(1+normal_slope ** 2)

        if edges_x == None or edges_y == None:
            print ("There is no critical region in this regime. Terminating calculation.")
            exit  ()

        phi_b_spin_min = 0.01
        phi_b_spin_max = 0.8

        phi_b_list = [np.linspace(phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1, mesh), \
                  np.linspace(phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3, mesh), \
                  np.linspace(phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5, mesh), \
                  np.linspace(phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9, mesh), \
                  np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh), \
                  np.logspace(-20, np.log10(0.9), mesh)] # , \
                  # np.hstack((np.linspace(phi_b_spin_min*0.001, phi_b_spin_min, mesh), np.linspace(phi_b_spin_max, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh)))] # \
                  # np.hstack((np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh), np.linspace(np.min(crits[:,1]), np.max(crits[:,1]), mesh)))]

        pool    = mp.Pool (processes=len(phi_b_list))
        results = pool.starmap(perform_sweep, zip(phi_b_list, itertools.repeat(mesh), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chi_ac), itertools.repeat(crits[0]), itertools.repeat(normal_to_crit), itertools.repeat(phi_a_edge), itertools.repeat(phi_b_edge) ) )
        pool.close ()
        pool.join  ()


        for idx,res in enumerate(results):
            if len(res) == 0:
                print (f"No critical condition in process {idx}.")
                continue
            for i in range(len (res[0]) ):
                f.write  ("{:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f} {:<1.20f}\n"\
        .format(res[2][i], res[0][i,0], res[0][i,1], 1-res[0][i,0]-res[0][i,1], res[1][i,0], res[1][i,1], 1-res[1][i,0]-res[1][i,1]))

    f.close ()

    stop = time.time ()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)


