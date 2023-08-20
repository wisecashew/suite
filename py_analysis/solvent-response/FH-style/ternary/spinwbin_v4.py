import numpy as np
import pandas as pd
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
from scipy.spatial import ConvexHull
import sys
import argparse
import time
import multiprocessing as mp
import os
import itertools
import warnings
import linecache

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

np.set_printoptions(threshold=sys.maxsize)
os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush()


import argparse
parser = argparse.ArgumentParser(description="Create a ternary /spin/odal and /bin/odal diagram.")
parser.add_argument('--chiac',        metavar='chi_ac', dest='chi_ac',       type=float,      action='store',             help='enter A-C exchange parameter.')
parser.add_argument('--chiab',        metavar='chi_ab', dest='chi_ab',       type=float,      action='store',             help='enter A-B exchange parameter.')
parser.add_argument('--chibc',        metavar='chi_bc', dest='chi_bc',       type=float,      action='store',             help='enter B-C exchange parameter.')
parser.add_argument('-N',             metavar='N',      dest='N',            type=int,        action='store',             help='degree of polymerization of B.')
parser.add_argument('--dumpfile',     dest='dumpfile',  type=str,            action='store',  help="name of file where the skeleton was dumped.")
parser.add_argument('--bin-boundary', dest='boundary',  type=str,            action='store',  help="name of file where you will dump out your solution for the binodal")
parser.add_argument('--tielines',     dest='tl',        action='store_true', default=False,   help="Option to include if you want to see all tie-lines.")
parser.add_argument('--image',        dest='img',       type=str,            action='store',  help="name of image generated.")
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
                            roots_down = np.vstack ((roots_down,r_tup))

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

################################

def sort_solutions (unsorted_upper, unsorted_lower, center, central_axis):

    top_half    = np.empty((0,3))
    bottom_half = np.empty((0,3))
    theta_1 = []
    theta_2 = []

    for s1, s2, in zip (unsorted_upper, unsorted_lower):
        d1 = (s1[0:2] - center)/np.linalg.norm(s1[0:2] - center)
        d2 = (s2[0:2] - center)/np.linalg.norm(s2[0:2] - center)
        clock1     = np.sign(np.cross (central_axis, d1))
        clock2     = np.sign(np.cross (central_axis, d2))

        if clock1 == 1:
            if clock2 == -1:
                pass
            else:
                print (f"This is strange. point1 = {s1}, point2 = {s2}.", flush=True)
                exit ()
            t1 = np.arccos (np.dot(d1, central_axis))
            theta_1.append (t1)
            top_half = np.vstack ((top_half, s1))
            t2 = np.arccos (np.dot(d2, central_axis))
            theta_2.append (t2)
            bottom_half = np.vstack ((bottom_half, s2))
        elif clock1 == -1:
            if clock2 == 1:
                pass
            else:
                print (f"This is strange. point1 = {s1}, point2 = {s2}.", flush=True)
                exit ()
            t1 = np.arccos (np.dot(d1, central_axis))
            theta_2.append (t1)
            bottom_half = np.vstack ((bottom_half, s1))
            t2 = np.arccos (np.dot(d2, central_axis))
            theta_1.append (t2)
            top_half    = np.vstack ((top_half, s2))

        elif clock1 == 0:
            if clock2 == 0:
                print (f"We are at crit point.", flush=True)
            else:
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
            t1          = np.arccos (np.dot(direction, central_axis))
            theta_1.append (t1)
            top_half    = np.vstack ((top_half, s1))
            t2          = np.arccos (np.dot(direction, central_axis))
            theta_2.append (t2)
            bottom_half = np.vstack ((bottom_half, s2))

        else:
            print ("Something's wrong.", flush=True)

    sorted_t1_idx = np.argsort (theta_1)
    top_half      = top_half [sorted_t1_idx]
    sorted_t2_idx = np.argsort (theta_2)
    bottom_half   = bottom_half[sorted_t2_idx]

    return top_half, bottom_half

###############################

def add_interpolated_rows (original_array, M):
    result_array = []

    for i in range(len(original_array) - 1):
        start_row = original_array[i]
        end_row = original_array[i + 1]
        step = (end_row - start_row) / (M + 1)

        interpolated_rows = [start_row + step * j for j in range(1, M + 1)]
        result_array.extend([start_row] + interpolated_rows)

    result_array.append(original_array[-1])
    return np.array(result_array)

#############################3

def add_rows_between_largest_gap (array, M):

    diff        = np.diff (array, axis=0)
    distances   = np.linalg.norm (diff, axis=1)
    max_dists   = np.argmax (distances)
    insert_rows = np.linspace (array[max_dists], array[max_dists+1], M)
    array       = np.insert (array, max_dists+1, insert_rows[1:-1],0)

    return array, max_dists

###############################

def max_dists_on_binodal (top_half, bottom_half):

    diff_top         = np.diff        (top_half, axis=0)
    distances        = np.linalg.norm (diff_top, axis=1)
    max_dists_top    = np.max         (distances)

    diff_bottom      = np.diff        (bottom_half, axis=0)
    distances        = np.linalg.norm (diff_bottom, axis=1)
    max_dists_bottom = np.max         (distances)

    max_dist = np.max ([max_dists_top, max_dists_bottom])
    
    return max_dist

######################################

def refined_binodal_v4 (side_1, side_2, nadded_rows):

    # together = np.vstack ((side_1, side_2))

    side_1, m1 = add_rows_between_largest_gap (side_1, nadded_rows)
    side_2, m2 = add_rows_between_largest_gap (side_2, nadded_rows)

    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v4...", flush=True)
    print (f"from {side_1[m1+1]} to {side_1[m1+nadded_rows-1]}", flush=True)
    # print (side_1[m1+1:m1+20])
    # print (side_2[m2+1:m2+20])

    for idx, pt in enumerate (side_1[m1+1:m1+nadded_rows-1]):
        if idx%25==0: print (f"idx = {idx} @ x,y = {pt[0],pt[1]}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        root_store = []
        dist_store = []

        for tidx, tpt in enumerate (side_2[m2+1+idx-nadded_rows:m2+1+idx+nadded_rows]):
            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]])

            # if the roots are "bad" roots, just write them out as bad
            if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
                continue

            else:
                fa = [pt[0], root[1]]
                fb = [root[0], root[2]]
                fc = [1-pt[0]-root[0], 1-root[1]-root[2]]
                p1 = np.array([pt[0], root[0], 1-root[0]-pt[0]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                continue

            elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                continue

            elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                continue

            # if the roots are basically the same point, write them out as bad 
            root_store.append ((p1,p2))
            dist_store.append (np.linalg.norm (p1-p2))

        # choose the root that is furthest away
        try:
            best_root = np.argmax (dist_store)
            root_combo = root_store[best_root]
            side_1[m1+1+idx] = root_combo[0]
            side_2[m1+1+idx] = root_combo[1]

        except:
            print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

    

    return [side_1, side_2]

######################################

def refined_binodal_v5 (side_1, side_2, nadded_rows):

    side_1, m1 = add_rows_between_largest_gap (side_1, nadded_rows)
    side_2, m2 = add_rows_between_largest_gap (side_2, nadded_rows)

    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v5...", flush=True)

    for idx, pt in enumerate (side_1[m1+1:m1+nadded_rows-1]):
        if idx%25==0: print (f"idx = {idx} @ x, y = {pt[0],pt[1]}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(phi[0], pt[1]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], pt[1]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], pt[1]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        root_store = []
        dist_store = []

        for tidx, tpt in enumerate (side_2[m2+1+idx-2*nadded_rows:m2+1+idx+2*nadded_rows]):
            root = fsolve (mu_equations, [pt[0], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if (np.abs(np.array(mu_equations(root))) > 1e-12).any():
                continue

            else:
                fa = [root[0], root[1]]
                fb = [pt[1]  , root[2]]
                fc = [1-root[0]-pt[1], 1-root[1]-root[2]]
                p1 = np.array([root[0], pt[1], 1-root[0]-pt[1]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                continue

            elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                continue

            elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                continue

            # if the roots are basically the same point, write them out as bad 
            root_store.append ((p1,p2))
            dist_store.append (np.linalg.norm (p1-p2))

        # choose the root that is furthest away
        try:
            best_root = np.argmax (dist_store)
            root_combo = root_store[best_root]
            side_1[m1+1+idx] = root_combo[0]
            side_2[m1+1+idx] = root_combo[1]

        except:
            print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

    return [side_1, side_2]

######################################

def root_finder_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string):

    if scale_string == "phi_a":
        scale_a = 1e+6
        scale_b = 1
    elif scale_string == "phi_b":
        scale_a = 1
        scale_b = 1e+6
    else:
        print (f"Bad string provided: {scale_string}.")
        exit()

    lower_guesses = np.linspace (sol_lower[max_ind], sol_lower[max_ind+1], 100)

    direction     = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_upper   = np.arccos(np.dot(direction, central_axis))
    theta1        = theta_upper [max_ind]
    theta2        = theta_upper [max_ind+1]

    to_keep       = (theta_upper > theta1) & (theta_upper < theta2)

    unsolved_upper_bin = binodal_upper[bad_idx, 0:2]
    unsolved_upper_bin = unsolved_upper_bin[to_keep]

    sol_bin_up    = np.empty((0,3))
    sol_bin_down  = np.empty((0,3))

    lower_guesses = np.linspace (sol_lower[max_ind], sol_lower[max_ind+1], 100)

    for idx in range (0, len(unsolved_upper_bin), 100):
        print (f"@ idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(phi[0], unsolved_upper_bin[idx][1]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], unsolved_upper_bin[idx][1]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], unsolved_upper_bin[idx][1]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        for iidx in range(len(lower_guesses)):
            root = fsolve (mu_equations, [unsolved_upper_bin[idx][0], lower_guesses[iidx][0]/scale_a, lower_guesses[iidx][1]/scale_b])

            if (np.abs(np.array(mu_equations(root)))>1e-6).any():
                continue

            else:
                p1 = np.array([root[0], unsolved_upper_bin[idx][1], 1-root[0]-unsolved_upper_bin[idx][1]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

                if np.linalg.norm(p1-p2) < 1e-6:
                    continue

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                    continue

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue
                else:
                    print ("HIT!", flush=True, end=' ')
                    print (f"p1 = {p1}, p2 = {p2}!", flush=True)
                    sol_bin_up   = np.vstack((sol_bin_up, p1))
                    sol_bin_down = np.vstack((sol_bin_down,p2))
                    break

    for idx in range (0, len(unsolved_upper_bin), 100):
        print (f"@ idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(unsolved_upper_bin[idx][0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(unsolved_upper_bin[idx][0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(unsolved_upper_bin[idx][0], phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        for iidx in range(len(lower_guesses)):
            root = fsolve (mu_equations, [unsolved_upper_bin[idx][1], lower_guesses[iidx][0]/scale_a, lower_guesses[iidx][1]/scale_b])

            if (np.abs(np.array(mu_equations(root)))>1e-6).any():
                continue

            else:
                p1 = np.array([unsolved_upper_bin[idx][0], root[0], 1-root[0]-unsolved_upper_bin[idx][0]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

                if np.linalg.norm(p1-p2) < 1e-6:
                    continue

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                    continue

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue
                else:
                    print ("HIT!", flush=True, end=' ')
                    print (f"p1 = {p1}, p2 = {p2}!", flush=True)
                    sol_bin_up   = np.vstack((sol_bin_up, p1))
                    sol_bin_down = np.vstack((sol_bin_down,p2))
                    break



    sol_upper              = np.vstack((sol_upper, sol_bin_up  ))
    sol_lower              = np.vstack((sol_lower, sol_bin_down))   

    # sort the solutions
    #########################
    direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_upper            = np.arccos (np.dot (direction, central_axis))
    sorted_theta_upper_idx = np.argsort (theta_upper)

    sol_upper = sol_upper [sorted_theta_upper_idx]
    sol_lower = sol_lower [sorted_theta_upper_idx]

    return sol_upper, sol_lower

######################################

def root_finder_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string):

    if scale_string == "phi_a":
        scale_a = 1e+6
        scale_b = 1
    elif scale_string == "phi_b":
        scale_a = 1
        scale_b = 1e+6
    else:
        print (f"Bad string provided: {scale_string}.")
        exit()

    upper_guesses = np.linspace (sol_upper[max_ind], sol_upper[max_ind+1], 100)

    direction     = (sol_lower[:,0:2] - center)/np.linalg.norm(sol_lower[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_lower   = np.arccos(np.dot(direction, central_axis))
    theta1        = theta_lower [max_ind]
    theta2        = theta_lower [max_ind+1]

    to_keep       = (theta_upper > theta1) & (theta_upper < theta2)

    unsolved_lower_bin = binodal_lower[bad_idx, 0:2]
    unsolved_lower_bin = unsolved_lower_bin[to_keep]

    sol_bin_up    = np.empty((0,3))
    sol_bin_down  = np.empty((0,3))

    upper_guesses = np.linspace (sol_upper[max_ind], sol_upper[max_ind+1], 100)

    for idx in range (0, len(unsolved_lower_bin), 100):
        print (f"@ idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(phi[0], unsolved_lower_bin[idx][1]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], unsolved_lower_bin[idx][1]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], unsolved_lower_bin[idx][1]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        for iidx in range(len(upper_guesses)):
            root = fsolve (mu_equations, [unsolved_lower_bin[idx][0], upper_guesses[iidx][0]/scale_a, upper_guesses[iidx][1]/scale_b])

            if (np.abs(np.array(mu_equations(root)))>1e-6).any():
                continue

            else:
                p1 = np.array([root[0], unsolved_lower_bin[idx][1], 1-root[0]-unsolved_lower_bin[idx][1]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

                if np.linalg.norm(p1-p2) < 1e-6:
                    continue

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                    continue

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue
                else:
                    print ("HIT!", flush=True, end=' ')
                    print (f"p1 = {p1}, p2 = {p2}!", flush=True)
                    sol_bin_up   = np.vstack((sol_bin_up, p1))
                    sol_bin_down = np.vstack((sol_bin_down,p2))
                    break

    for idx in range (0, len(unsolved_lower_bin), 100):
        print (f"@ idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(unsolved_lower_bin[idx][0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(unsolved_lower_bin[idx][0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(unsolved_lower_bin[idx][0], phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        for iidx in range(len(upper_guesses)):
            root = fsolve (mu_equations, [unsolved_lower_bin[idx][1], upper_guesses[iidx][0]/scale_a, upper_guesses[iidx][1]/scale_b])

            if (np.abs(np.array(mu_equations(root)))>1e-6).any():
                continue

            else:
                p1 = np.array([unsolved_lower_bin[idx][0], root[0], 1-root[0]-unsolved_lower_bin[idx][0]])
                p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

                if np.linalg.norm(p1-p2) < 1e-6:
                    continue

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                    continue

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                    continue
                else:
                    print ("HIT!", flush=True, end=' ')
                    print (f"p1 = {p1}, p2 = {p2}!", flush=True)
                    sol_bin_up   = np.vstack((sol_bin_up, p1))
                    sol_bin_down = np.vstack((sol_bin_down,p2))
                    break



    sol_upper              = np.vstack((sol_upper, sol_bin_up  ))
    sol_lower              = np.vstack((sol_lower, sol_bin_down))

    # sort the solutions
    #########################
    direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_upper            = np.arccos (np.dot (direction, central_axis))
    sorted_theta_upper_idx = np.argsort (theta_upper)

    sol_upper = sol_upper [sorted_theta_upper_idx]
    sol_lower = sol_lower [sorted_theta_upper_idx]

    return sol_upper, sol_lower

######################################

def binodal_plotter (fig, ax, dumpfile, chi_ab, chi_bc, chi_ac, va, vb, vc, crit_points):

    try:
        df = pd.read_csv (dumpfile, sep='\s+', engine="python", skiprows=1, names=["dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2"])
        df = df.loc[df["dmu"]<5]
    except FileNotFoundError:
        print (f"File called {dumpfile} was not found. This was likely because pskelbin.py could not find reasonable guesses. Please check your parameters and inputs and try again.", flush=True)
        exit ()

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    # original point and the guess for the root... 
    phi_a_upper = df["phi_a1"].values; phi_a_lower = df["phi_a2"].values
    phi_b_upper = df["phi_b1"].values; phi_b_lower = df["phi_b2"].values
    phi_c_upper = df["phi_c1"].values; phi_c_lower = df["phi_c2"].values

    binodal_upper  = np.zeros ((phi_a_upper.shape[0],3))
    binodal_lower  = np.zeros ((phi_a_upper.shape[0],3))

    sol_upper = np.empty((0,3))
    sol_lower = np.empty((0,3))

    bad_idx  = []
    good_idx = []
    # f = open (args.boundary, 'w')
    print ("Start processing the dumpfile and find roots.", flush=True)
    print (f"I will be processing an array of size = {phi_a_upper.size}.", flush=True)

    for idx in range (len(phi_a_upper)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b_upper[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b_upper[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b_upper[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a_upper[idx], phi_a_lower[idx], phi_b_lower[idx]])
        binodal_upper [idx, :] = np.array([phi_a_upper[idx], phi_b_upper[idx], 1-phi_a_upper[idx]-phi_b_upper[idx]])
        binodal_lower [idx, :] = np.array([phi_a_lower[idx], phi_b_lower[idx], 1-phi_a_lower[idx]-phi_b_lower[idx]])

        # if the roots are "bad" roots, just write them out as bad
        if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
            bad_idx.append (idx)
            continue

        else:
            fa = [root[0], root[1]]
            fb = [phi_b_upper[idx], root[2]]
            fc = [1-root[0]-phi_b_upper[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b_upper[idx], 1-root[0]-phi_b_upper[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if np.linalg.norm(p1-p2) < 1e-6:
                bad_idx.append (idx)

            elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                bad_idx.append (idx)

            elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                bad_idx.append (idx)

            elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                bad_idx.append (idx)

            else:
                good_idx.append (idx)
                sol_upper = np.vstack ((sol_upper,p1))
                sol_lower = np.vstack ((sol_lower,p2))

    # start partitioning along a certain axis
    center         = np.mean (crit_points, axis=0)[:2]
    central_axis   = (crit_points[0,:2]-center)/np.linalg.norm (crit_points[0,:2]-center)

    print (f"center = {center}", flush=True)
    print (f"central axis = {central_axis}", flush=True)

    direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_upper            = np.arccos (np.dot (direction, central_axis))
    sorted_theta_upper_idx = np.argsort (theta_upper)
    theta_upper            = theta_upper[sorted_theta_upper_idx]

    direction              = (sol_lower[:,0:2] - center)/np.linalg.norm(sol_lower[:,0:2] - center, axis=1)[:, np.newaxis]
    theta_lower            = np.arccos (np.dot (direction, central_axis))
    sorted_theta_lower_idx = np.argsort (theta_lower)
    theta_lower            = theta_lower[sorted_theta_lower_idx]

    sol_upper              = sol_upper[sorted_theta_upper_idx]
    sol_lower              = sol_lower[sorted_theta_lower_idx]

    # I have sorted the solutions 
    # now, if the solution curve is sufficiently close, no need to perform more detailed searches -- so find maximum distances between points on the solution curves
    # find differences between lower and upper curves

    diff_up       = np.linalg.norm(sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
    max_diff_up   = np.max(diff_up)

    diff_down     = np.linalg.norm(sol_lower[1:][:,0:2] - sol_lower[:-1][:,0:2], axis=1)
    max_diff_down = np.max(diff_down)

    while (max_diff_up > 0.1 or max_diff_down > 0.1):
        
        # if max_diff_up > max_diff_down:
        print (f"Running a search on the top half...", flush=True)
        # start running a finer search 
        difference    = np.linalg.norm (sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
        max_ind       = np.argmax (difference)

        #upper edges 
        u1 = sol_upper[max_ind][0:2]
        u2 = sol_upper[max_ind+1][0:2]
        umin = np.min([u1,u2])
        l1 = sol_lower[max_ind][0:2]
        l2 = sol_lower[max_ind+1][0:2]
        lmin = np.min([l1,l2])

        if umin > lmin:
            # scale lower guesses
            if   ( abs(l1[0] - l2[0]) < abs(l1[1]-l2[1]) ):
                scale_string = "phi_a"
                sol_upper, sol_lower = root_finer_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string)
        

            elif ( abs(l1[0] - l2[0]) > abs(l1[1]-l2[1]) ):
                scale_string = "phi_b"
                sol_upper, sol_lower = root_finer_with_scaling_lower (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string)

        elif umin < lmin:
            # scale upper guesses
            if   ( abs(u1[0] - u2[0]) < abs(u1[1]-u2[1]) ):
                scale_string = "phi_a"
                sol_upper, sol_lower = root_finer_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string)

            elif ( abs(u1[0] - u2[0]) > abs(u1[1]-u2[1]) ):
                scale_string = "phi_b"
                sol_upper, sol_lower = root_finer_with_scaling_upper (sol_upper, sol_lower, max_ind, binodal_upper, binodal_lower, bad_idx, center, central_axis, scale_string)


        diff_up       = np.linalg.norm(sol_upper[1:] - sol_upper[:-1], axis=1)
        max_diff_up   = np.max(diff_up)

        diff_down     = np.linalg.norm(sol_lower[1:] - sol_lower[:-1], axis=1)
        max_diff_down = np.max(diff_down)


    ##################################################

    print ("Broke out of initial solver. Time to refine...", flush=True)
    # WE NOW HAVE A PRETTY CLEANED OUT BINODAL. 
    # ALL THAT REMAINS IS TO REFINE IT TO MAKE IT LOOK NICE.

    diff_up       = np.linalg.norm(sol_upper[1:][:,0:2] - sol_upper[:-1][:,0:2], axis=1)
    max_diff_up   = np.max(diff_up)

    diff_down     = np.linalg.norm(sol_lower[1:][:,0:2] - sol_lower[:-1][:,0:2], axis=1)
    max_diff_down = np.max(diff_down)

    print ("Being refining!", flush=True)
    while (max_diff_up > 0.01 or max_diff_down > 0.01):
        print (f"@ max_diff_up = {max_diff_up}, max_diff_down = {max_diff_down}...")
        nadded_rows = 100
        sol_upper, sol_lower = refined_binodal_v4 (sol_upper, sol_lower, nadded_rows)

        direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
        theta_upper            = np.arccos (np.dot (direction, central_axis))
        sorted_theta_upper_idx = np.argsort (theta_upper)
        theta_upper            = theta_upper[sorted_theta_upper_idx]
        sol_upper              = sol_upper[sorted_theta_upper_idx]
        sol_lower              = sol_lower[sorted_theta_upper_idx]

        sol_upper, sol_lower = refined_binodal_v5 (sol_upper, sol_lower, nadded_rows)
        direction              = (sol_upper[:,0:2] - center)/np.linalg.norm(sol_upper[:,0:2] - center, axis=1)[:, np.newaxis]
        theta_upper            = np.arccos (np.dot (direction, central_axis))
        sorted_theta_upper_idx = np.argsort (theta_upper)
        theta_upper            = theta_upper[sorted_theta_upper_idx]
        sol_upper              = sol_upper[sorted_theta_upper_idx]
        sol_lower              = sol_lower[sorted_theta_upper_idx]

        diff_up       = np.linalg.norm(sol_upper[1:] - sol_upper[:-1], axis=1)
        max_diff_up   = np.max(diff_up)

        diff_down     = np.linalg.norm(sol_lower[1:] - sol_lower[:-1], axis=1)
        max_diff_down = np.max(diff_down)


    ref_bin = [sol_upper, sol_lower]
    print ("Both crit points should be well-populated.", flush=True)

    # this is the binodal
    ax.scatter (ref_bin[0][:,0], ref_bin[0][:,1], c='k', s=0.125, zorder=11)
    ax.scatter (ref_bin[1][:,0], ref_bin[1][:,1], c='k', s=0.125, zorder=11)

    if args.tl:
        for i in range (len(ref_bin[0])):
            ax.plot    ([ref_bin[0][i,0],ref_bin[1][i,0]], [ref_bin[0][i,1],ref_bin[1][i,1]], lw=0.5, ls='--', markersize=0, zorder=10)

    ff = open (args.boundary, 'w')
    ff.write ("phi_s_top|phi_p_top|phi_c_top|phi_s_bot|phi_p_bot|phi_c_bot\n")
    ff.write (f"{crit_points[0,0]}|{crit_points[0,1]}|{1-crit_points[0,0]-crit_points[0,1]}|{crit_points[1,0]}|{crit_points[1,1]}|{1-crit_points[1,0]-crit_points[1,1]}\n")
    for i in range (len(ref_bin[0])):
        ff.write (f"{ref_bin[0][i][0]}|{ref_bin[0][i][1]}|{ref_bin[0][i][2]}|{ref_bin[1][i][0]}|{ref_bin[1][i][1]}|{ref_bin[1][i][2]}\n")

    return

############################


if __name__=="__main__":

    start = time.time()

    ###########################
    N = args.N
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc
    chi_ac = args.chi_ac
    dumpfile = args.dumpfile
    ############################

    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(3,3))
    ax  = plt.axes ()

    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4* N * (1 - 2* phi_s * chi_sc + 2 * phi_s ** 2 * chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    (-1 + 2 * phi_s * chi_sc + N * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc:  1 / (2*N * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: 1 - 2 * phi_s * chi_sc + N * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )


    roots_up, roots_down = find_crit_point (N, chi_ac, chi_ab, chi_bc)
    ax.scatter (roots_up[:,0]  , roots_up[:,1]  , color='k', edgecolors='steelblue', s=0.1, zorder=11)
    ax.scatter (roots_down[:,0], roots_down[:,1], color='k', edgecolors='steelblue', s=0.1, zorder=11)

    crits = np.vstack ((roots_up, roots_down))
    if len(crits) == 2:
        pass
    else:
        print ("Number of critical points = ",len(crits), flush=True)
    print ("critical points = \n",crits, flush=True)

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    print ("Begin creating the meshes and painting the ternary diagram...", flush=True)
    p_a_space = np.arange (0.001, 1-0.001, 0.001)
    p_a = np.repeat (p_a_space, len(p_a_space))

    p_b = np.zeros (p_a.shape)
    for i in range (len(p_a_space)):
        p_b[i*len(p_a_space):(i+1)*len(p_a_space)] = np.linspace (0.001, 1-p_a_space[i], len(p_a_space))

    vals = stab_crit (p_a, p_b, chi_ab, chi_bc, chi_ac)

    to_keep = ~np.isnan(vals)

    vals = vals[to_keep]
    p_a  = p_a [to_keep]
    p_b  = p_b [to_keep]

    vmax = np.max (vals)
    vmin = np.min (vals)

    if np.sign (vmax) == np.sign (vmin):
        if np.sign (vmax) >=0:
            vmin = -vmax
        else:
            vmax = -vmin

    print (f"vmin = {vmin}, vmax = {vmax}", flush=True)

    norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
    cols = cm.bwr (norm (vals))

    # Plot the points
    p_c = 1 - p_a - p_b
    ax.scatter  (p_a, p_b, s=0.01, color=cols, zorder=0, clip_on=True)
    print ("Painted the ternary diagram!", flush=True)


    ax.set_xlabel('$\\phi _{S}$')
    ax.set_ylabel('$\\phi _{P}$')
    print ("We have plotted the spinodal region!\n", flush=True)
    print ("###########################################################\n", flush=True)
    print ("Start binodal plotting...\n", flush=True)

    va = 1
    vb = N
    vc = 1

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    binodal_plotter (fig, ax, dumpfile, chi_ab, chi_bc, chi_ac, va, vb, vc, crits)
    print ("Done with binodal plotting!", flush=True)

    ax.plot (crits[:,0], crits[:,1], c='darkred', zorder=12, lw=0.5)

    # Set axis limits
    ax.set_xlim (-0.01, 1)
    ax.set_ylim (-0.01, 1)

    ax.grid()

    fig.savefig  (args.img, dpi=1200)
    print ("Completed heat map computation.", flush=True)
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.", flush=True)
