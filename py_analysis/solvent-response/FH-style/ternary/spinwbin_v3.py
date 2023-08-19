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
parser.add_argument('--nproc',        dest='nproc',     type=int,            action='store',  help="number of processes to gather.")
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
        root = fsolve (send_to_fsolve_r1, g, xtol=1e-12)

        if abs(send_to_fsolve_r1(root)) < 1e-12:

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
        root = fsolve (send_to_fsolve_r2, g, xtol=1e-12)

        if abs(send_to_fsolve_r2(root)) < 1e-12:

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
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
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
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
                exit ()
            t1 = np.arccos (np.dot(d1, central_axis))
            theta_2.append (t1)
            bottom_half = np.vstack ((bottom_half, s1))
            t2 = np.arccos (np.dot(d2, central_axis))
            theta_1.append (t2)
            top_half    = np.vstack ((top_half, s2))

        elif clock1 == 0:
            if clock2 == 0:
                print (f"We are at crit point.")
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

###############################

def refined_binodal (side_1, side_2):

    together = np.vstack ((side_1, side_2))
    bin1 = np.empty ((0,3))
    bin2 = np.empty ((0,3))
    print (f"together = {together.shape}", flush=True)
    for idx, pt in enumerate (together):
        print (f"idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]
        root_store = []
        dist_store = []
        for tidx, tpt in enumerate (together):
            if idx == tidx:
                continue

            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
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
        best_root = np.argmax (dist_store)
        root_combo = root_store[best_root]
        bin1 = np.vstack((bin1, root_combo[0]))
        bin2 = np.vstack((bin2, root_combo[1]))

    return [bin1, bin2]

######################################

def refined_binodal_v2 (side_1, side_2):

    # together = np.vstack ((side_1, side_2))
    bin1 = np.empty ((0,3))
    bin2 = np.empty ((0,3))
    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal...", flush=True)
    for idx, pt in enumerate (side_1):
        print (f"idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]
        root_store = []
        dist_store = []
        for tidx, tpt in enumerate (side_2):

            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
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
            bin1 = np.vstack((bin1, root_combo[0]))
            bin2 = np.vstack((bin2, root_combo[1]))
        except:
            print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)

    for idx, pt in enumerate (side_2):
        print (f"idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]
        root_store = []
        dist_store = []
        for tidx, tpt in enumerate (side_1):

            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
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
            bin1 = np.vstack((bin1, root_combo[0]))
            bin2 = np.vstack((bin2, root_combo[1]))
        except:
            print (f"Problem with a particular value. No solution found for pt = {pt}", flush=True)


    return [bin1, bin2]

######################################

def refined_binodal_v3 (side_1, side_2, added_rows):

    # together = np.vstack ((side_1, side_2))
    bin1 = np.empty ((0,3))
    bin2 = np.empty ((0,3))
    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal...", flush=True)
    for idx, pt in enumerate (side_1):
        print (f"idx = {idx}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        root_store = []
        dist_store = []
        counter = idx // added_rows 
        if counter < 6:
             lower = 0
             upper = counter + 10
        else:
             lower = counter - 5
             upper = counter + 5
        for tidx, tpt in enumerate (side_2[(added_rows)*(lower):(added_rows)*(upper)]):

            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if (np.abs(np.array(mu_equations(root))) > 1e-12).any():
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
            best_root  = np.argmax (dist_store)
            root_combo = root_store[best_root]
            bin1 = np.vstack((bin1, root_combo[0]))
            bin2 = np.vstack((bin2, root_combo[1]))
        except:
            print (f"Problem with a particular value. No solution found for pt = {pt}")

    return [bin1, bin2]

######################################

def refined_binodal_v4 (side_1, side_2, added_rows):

    # together = np.vstack ((side_1, side_2))

    side_1, m1 = add_rows_between_largest_gap (side_1, added_rows)
    side_2, m2 = add_rows_between_largest_gap (side_2, added_rows)

    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v4...", flush=True)
    print (f"from {side_1[m1+1]} to {side_1[m1+added_rows-1]}", flush=True)
    # print (side_1[m1+1:m1+20])
    # print (side_2[m2+1:m2+20])

    for idx, pt in enumerate (side_1[m1+1:m1+added_rows-1]):
        if idx%25==0: print (f"idx = {idx} @ x,y = {pt[0],pt[1]}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(pt[0], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(pt[0], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(pt[0], phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        root_store = []
        dist_store = []

        for tidx, tpt in enumerate (side_2[m2+1+idx-added_rows:m2+1+idx+added_rows]):
            root = fsolve (mu_equations, [pt[1], tpt[0], tpt[1]], xtol=1e-12)

            # if the roots are "bad" roots, just write them out as bad
            if (np.abs(np.array(mu_equations(root))) > 1e-12).any():
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

def refined_binodal_v5 (side_1, side_2, added_rows):

    # together = np.vstack ((side_1, side_2))

    side_1, m1 = add_rows_between_largest_gap (side_1, added_rows)
    side_2, m2 = add_rows_between_largest_gap (side_2, added_rows)

    print (f"side_1.shape = {side_1.shape}, side_2.shape = {side_2.shape}.\nRefining binodal with v5...", flush=True)
    # print (f"from {side_1[m1+1]} to {side_1[m1+added_rows]}")
    # print (side_1[m1+1:m1+20])
    # print (side_2[m2+1:m2+20])

    for idx, pt in enumerate (side_1[m1+1:m1+added_rows-1]):
        if idx%25==0: print (f"idx = {idx} @ x, y = {pt[0],pt[1]}...", flush=True)
        def mu_equations (phi):
            eq1 = mu_a(phi[0], pt[1]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], pt[1]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], pt[1]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        root_store = []
        dist_store = []

        for tidx, tpt in enumerate (side_2[m2+1+idx-2*added_rows:m2+1+idx+2*added_rows]):
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

# go_through_indices goes through every single "bad" solution I have, 
# and tries to find a good solution for it using guesses created
# that are actual points on the binodal.

def go_through_indices (bad_idx_subset, good_idx, initg, bincloser, binfurther):
    sol1_bg = np.empty ((0,3))
    sol2_bg = np.empty ((0,3))
    for i, bidx in enumerate (bad_idx_subset):
        phi_a  = initg[bidx, 0]
        def mu_equations (phi):
            eq1 = mu_a(phi_a, phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi_a, phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi_a, phi[0]) - mu_c(phi[1], phi[2])
            return [eq1, eq2, eq3]

        # go through good indices 
        for j, gidx in enumerate (good_idx):
            root = fsolve (mu_equations, [bincloser[gidx, 1], binfurther[gidx, 0], binfurther[gidx,1]], xtol=1e-12)
            p1   = np.array ([phi_a, root[0], 1-phi_a-root[0]])
            p2   = np.array ([root[1], root[2], 1-root[1]-root[2]])
            if (np.abs ( np.array (mu_equations (root))) > 1e-12).any ():
                continue
            elif (np.linalg.norm(sol1_bg - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2_bg - p2, axis=-1)<1e-6).any():
                continue

            elif np.linalg.norm (p1-p2) < 1e-6:
                continue

            else:
                sol1_bg = np.vstack ((sol1_bg, p1))
                sol2_bg = np.vstack ((sol2_bg, p2))
                break

    return [sol1_bg, sol2_bg]


def binodal_plotter (fig, ax, dumpfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc, crit_points):

    try:
        df = pd.read_csv (dumpfile, sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])
        df = df.loc[df["dmu"]<1]
    except FileNotFoundError:
        print (f"File called {dumpfile} was not found. This was likely because pskelbin.py could not find reasonable guesses. Please check your parameters and inputs and try again.", flush=True)
        exit ()

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    # original point and the guess for the root... 
    phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
    phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
    phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

    init         = np.array ([phi_a,phi_b,phi_c]).T
    seed         = np.array ([phi_an, phi_bn, phi_cn]).T
    binodal_closer  = np.zeros ((phi_a.shape[0],3))
    binodal_further = np.zeros ((phi_a.shape[0],3))


    sol1 = np.empty((0,3))
    sol2 = np.empty((0,3))
    bad_idx  = []
    good_idx = []
    # f = open (args.boundary, 'w')
    print ("Start processing the dumpfile and find roots.", flush=True)
    for idx in range (len(phi_a)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]], xtol=1e-12)
        binodal_closer [idx, :] = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
        binodal_further[idx, :] = np.array([root[1], root[2], 1-root[1]-root[2]])
        # print (f"phi_b[{idx}] = {phi_b[idx]}")

        # if the roots are "bad" roots, just write them out as bad
        if ( np.abs(np.array(mu_equations(root))) > 1e-12).any():
            bad_idx.append (idx)
        else:
            fa = [root[0], root[1]]
            fb = [phi_b[idx], root[2]]
            fc = [1-root[0]-phi_b[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            # if the roots are basically the same point, write them out as bad 
            if np.linalg.norm (p1-p2) < 1e-6:
                continue

            else:
                # if one of the good solutions is same as something before, write it as bad 
                if (np.linalg.norm(sol1 - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2 - p2, axis=-1)<1e-6).any():
                    bad_idx.append (idx)

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                     bad_idx.append (idx)

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)

                else:
                    good_idx.append (idx)
                    sol1 = np.vstack ((sol1,p1))
                    sol2 = np.vstack ((sol2,p2))

    print   ("Everything has been written out. Start doing a bigger processing - with parallelization.", flush=True)
    # now, time to slowly convert the "bad" points to "good" points, in a parallel fashion
    pool  = mp.Pool ( processes=nproc )
    sol1_bg = np.empty ((0,3))
    sol2_bg = np.empty ((0,3))
    print (f"Spawning {nproc} processes...", flush=True)
    # divide up bad_idx 
    bad_idx_subsets = [bad_idx[i:i+nproc] for i in range (0, len(bad_idx), nproc)]
    print (f"About to send off processes...", flush=True)
    results = pool.starmap (go_through_indices, zip (bad_idx_subsets, itertools.repeat (good_idx), itertools.repeat (init), itertools.repeat (binodal_closer), itertools.repeat (binodal_further)))
    print (f"Processes completed! compiling results...", flush=True)

    for sols in results:
        sol1_bg = np.vstack ((sol1_bg, sols[0]))
        sol2_bg = np.vstack ((sol2_bg, sols[1]))

    pool.close ()
    pool.join  ()

    # FILTER NUMBER 1: GET RID OF ROWS THAT ARE TOO SIMILAR

    # curate sol1 and sol2
    sol1_bg, inds = np.unique (sol1_bg, axis=0, return_index=True)
    sol2_bg       = sol2_bg [inds, :]
    sol1_bg       = np.vstack((sol1_bg,sol1))
    sol2_bg       = np.vstack((sol2_bg,sol2))

    # collect all the points in sol1_bg and sol2_bg, and sort them for further processing
    sol_net = np.vstack ((sol1_bg, sol2_bg))

    # start partitioning along a certain axis
    center         = np.mean (crit_points, axis=0)[:2]
    central_axis   = (crit_points[0,:2]-center)/np.linalg.norm (crit_points[0,:2]-center)

    print (f"center = {center}", flush=True)
    print (f"central axis = {central_axis}", flush=True)

    side_1  = np.empty ((0,3))
    side_2  = np.empty ((0,3))
    theta_1 = []
    theta_2 = []

    for pts in sol_net:
        direction = (pts[0:2] - center)/np.linalg.norm(pts[0:2] - center)
        clock     = np.sign(np.cross (central_axis, direction))

        if clock == 1:
            t1 = np.arccos (np.dot(direction, central_axis))
            theta_1.append (t1)
            side_1 = np.vstack ((side_1, pts))
        elif clock == -1:
            t2 = np.arccos (np.dot(direction, central_axis))
            theta_2.append (t2)
            side_2 = np.vstack ((side_2, pts))

        elif clock == 0:
            t1 = np.arccos (np.dot(direction, central_axis))
            theta_1.append (t1)
            side_1 = np.vstack ((side_1, pts))

        else:
            print ("Something's wrong.", flush=True)

    theta_1 = np.array(theta_1)
    theta_2 = np.array(theta_2)
    sorted_t1_idx = np.argsort (theta_1)
    sorted_t2_idx = np.argsort (theta_2)

    print ("Roots have been found and sorted.", flush=True)

    # WE HAVE NOW SPLIT UP THE UNSTABLE REGION INTO TWO ZONES BY USING THE CRITICAL POINT
    # NOW THAT WE HAVE TWO ZONES AS TWO HALVES:side_1 and side_2. LETS USE THEM TO FIND THE BINODAL

    # now that we have divided into two sectors, let's sort them out 
    side_1 = side_1 [sorted_t1_idx]
    side_2 = side_2 [sorted_t2_idx]

    # I now have two sides of the binodal, somewhat arbitrarily split
    # now it is time to fill in the gaps
    if len(side_1) < 100:
        side_1 = add_interpolated_rows (side_1, 5)
    if len(side_2) < 100:
        side_2 = add_interpolated_rows (side_1, 5)
    print (f"About to start refining the binodal with v2...", flush=True)
    ref_bin = refined_binodal_v2 (side_1, side_2)
    diff  = ref_bin[0] - ref_bin[1]
    norms = np.linalg.norm (diff, axis=1)
    valid = np.where (norms >= 1e-6)[0]
    ref_bin[0] = ref_bin[0][valid,:]
    ref_bin[1] = ref_bin[1][valid,:]


    # I have used fsolve to find the first set of solution points. This is going to be rather sparse. 
    # what we do now is to fill up the empty spaces using finer searches
    top_half    = np.empty ((0,3))
    bottom_half = np.empty ((0,3))
    theta_1 = []
    theta_2 = []

    for s1, s2 in zip (ref_bin[0], ref_bin[1]):
        d1 = (s1[0:2] - center)/np.linalg.norm(s1[0:2] - center)
        d2 = (s2[0:2] - center)/np.linalg.norm(s2[0:2] - center)
        clock1     = np.sign(np.cross (central_axis, d1))
        clock2     = np.sign(np.cross (central_axis, d2))

        if clock1 == 1:
            if clock2 == -1:
                pass
            else:
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
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
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
                exit ()
            t1 = np.arccos (np.dot(d1, central_axis))
            theta_2.append (t1)
            bottom_half = np.vstack ((bottom_half, s1))
            t2 = np.arccos (np.dot(d2, central_axis))
            theta_1.append (t2)
            top_half    = np.vstack ((top_half, s2))

        elif clock1 == 0:
            if clock2 == 0:
                print (f"We are at crit point.")
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

    # print (f"theta = \n{theta_1}")
    sorted_t1_idx = np.argsort (theta_1)
    top_half      = top_half [sorted_t1_idx]
    sorted_t2_idx = np.argsort (theta_2)
    bottom_half   = bottom_half[sorted_t2_idx]

    # We have once again sorted the two solutions. Time to perform another refinement.
    # we will hunt for points along the y-axis

    print ("Beginning the v3 refinement...", flush=True)
    top_half    = add_interpolated_rows (top_half,    10)
    bottom_half = add_interpolated_rows (bottom_half, 10)

    ref_bin = refined_binodal_v3 (top_half, bottom_half, 10)
    diff  = ref_bin[0] - ref_bin[1]
    norms = np.linalg.norm (diff, axis=1)
    valid = np.where (norms >= 1e-6)[0]
    ref_bin[0] = ref_bin[0][valid,:]
    ref_bin[1] = ref_bin[1][valid,:]
    print (f"Refined binodal!", flush=True)

    top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)

    max_dist = max_dists_on_binodal (top_half, bottom_half)
    print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)

    refinement_iter = -1
    while max_dist > 0.1:
        # we will run another refinement.
        # we will hunt for points along the x-axis
        print ("Begin the v4 refinement...", flush=True)
        ref_bin     = refined_binodal_v4 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)


        # running another final refinement in the region that requires most space-filling
        # we will hunt for points along the y-axis
        print ("Begin the v5 refinement...", flush=True)
        ref_bin     = refined_binodal_v5 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)

        max_dist = max_dists_on_binodal (top_half, bottom_half)
        refinement_iter += 1
        print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)
        print (f"Refinement iteration: {refinement_iter}.", flush=True)

    # add the critical point to top_half and bottom_half
    print ("The binodal has been sufficiently filled up. Time to stack things around the critical points.", flush=True)
    crit_top_most = np.array([crit_points[0][0], crit_points[0][1], 1-crit_points[0][0]-crit_points[0][1]])
    crit_bot_most = np.array([crit_points[1][0], crit_points[1][1], 1-crit_points[1][0]-crit_points[1][1]])
    print (f"crit_top_most = {crit_top_most}", flush=True)
    print (f"crit_bot_most = {crit_bot_most}", flush=True)
    
    top_half    = np.vstack ((crit_top_most,top_half))
    bottom_half = np.vstack ((crit_top_most,bottom_half))

    max_dist = max_dists_on_binodal (top_half, bottom_half)
    print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)

    refinement_iter = -1
    while max_dist > 0.01:
        # we will run another refinement.
        # we will hunt for points along the x-axis
        print ("Begin the v4 refinement...", flush=True)
        ref_bin     = refined_binodal_v4 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)


        # running another final refinement in the region that requires most space-filling
        # we will hunt for points along the y-axis
        print ("Begin the v5 refinement...", flush=True)
        ref_bin     = refined_binodal_v5 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)

        max_dist = max_dists_on_binodal (top_half, bottom_half)
        refinement_iter += 1
        if refinement_iter == 4:
            break
        print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)
        print (f"Refinement iteration: {refinement_iter}.", flush=True)


    print ("One crit point should be well-populated.", flush=True)
    top_half    = np.vstack ((top_half   , crit_bot_most))
    bottom_half = np.vstack ((bottom_half, crit_bot_most))

    max_dist = max_dists_on_binodal (top_half, bottom_half)
    print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)

    refinement_iter = -1
    while max_dist > 0.01:
        # we will run another refinement.
        # we will hunt for points along the x-axis
        print ("Begin the v4 refinement...", flush=True)
        ref_bin     = refined_binodal_v4 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)


        # running another final refinement in the region that requires most space-filling
        # we will hunt for points along the y-axis
        print ("Begin the v5 refinement...", flush=True)
        ref_bin     = refined_binodal_v5 (top_half, bottom_half, 100)
        diff        = ref_bin[0] - ref_bin[1]
        norms       = np.linalg.norm (diff, axis=1)
        valid       = np.where (norms >= 1e-6)[0]
        ref_bin[0]  = ref_bin[0][valid,:]
        ref_bin[1]  = ref_bin[1][valid,:]
        top_half, bottom_half = sort_solutions (ref_bin[0], ref_bin[1], center, central_axis)

        max_dist = max_dists_on_binodal (top_half, bottom_half)
        refinement_iter += 1
        if refinement_iter == 4:
            break
        print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)
        print (f"Refinement iteration: {refinement_iter}.", flush=True)

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

def binodal_plotter_prototype (fig, ax, dumpfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc, crit_points):

    try:
        df = pd.read_csv (dumpfile, sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])
        df = df.loc[df["dmu"]<1]
    except FileNotFoundError:
        print (f"File called {dumpfile} was not found. This was likely because pskelbin.py could not find reasonable guesses. Please check your parameters and inputs and try again.", flush=True)
        exit ()

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    # original point and the guess for the root... 
    phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
    phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
    phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

    init         = np.array ([phi_a,phi_b,phi_c]).T
    seed         = np.array ([phi_an, phi_bn, phi_cn]).T
    binodal_closer  = np.zeros ((phi_a.shape[0],3))
    binodal_further = np.zeros ((phi_a.shape[0],3))


    sol1 = np.empty((0,3))
    sol2 = np.empty((0,3))
    bad_idx  = []
    good_idx = []
    # f = open (args.boundary, 'w')
    print ("Start processing the dumpfile and find roots.", flush=True)
    for idx in range (len(phi_a)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]], xtol=1e-12)
        binodal_closer [idx, :] = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
        binodal_further[idx, :] = np.array([root[1], root[2], 1-root[1]-root[2]])
        # print (f"phi_b[{idx}] = {phi_b[idx]}")

        # if the roots are "bad" roots, just write them out as bad
        if ( np.abs(np.array(mu_equations(root))) > 1e-12).any():
            bad_idx.append (idx)
        else:
            fa = [root[0], root[1]]
            fb = [phi_b[idx], root[2]]
            fc = [1-root[0]-phi_b[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            # if the roots are basically the same point, write them out as bad 
            if np.linalg.norm (p1-p2) < 1e-6:
                continue

            else:
                # if one of the good solutions is same as something before, write it as bad 
                if (np.linalg.norm(sol1 - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2 - p2, axis=-1)<1e-6).any():
                    bad_idx.append (idx)

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                     bad_idx.append (idx)

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)

                else:
                    good_idx.append (idx)
                    sol1 = np.vstack ((sol1,p1))
                    sol2 = np.vstack ((sol2,p2))

    print   ("Everything has been written out. Start doing a bigger processing - with parallelization.", flush=True)
    # now, time to slowly convert the "bad" points to "good" points, in a parallel fashion
    pool  = mp.Pool ( processes=nproc )
    sol1_bg = np.empty ((0,3))
    sol2_bg = np.empty ((0,3))
    print (f"Spawning {nproc} processes...", flush=True)
    # divide up bad_idx 
    bad_idx_subsets = [bad_idx[i:i+nproc] for i in range (0, len(bad_idx), nproc)]
    print (f"About to send off processes...", flush=True)
    results = pool.starmap (go_through_indices, zip (bad_idx_subsets, itertools.repeat (good_idx), itertools.repeat (init), itertools.repeat (binodal_closer), itertools.repeat (binodal_further)))
    print (f"Processes completed! compiling results...", flush=True)

    for sols in results:
        sol1_bg = np.vstack ((sol1_bg, sols[0]))
        sol2_bg = np.vstack ((sol2_bg, sols[1]))

    pool.close ()
    pool.join  ()

    # FILTER NUMBER 1: GET RID OF ROWS THAT ARE TOO SIMILAR

    # curate sol1 and sol2
    sol1_bg, inds = np.unique (sol1_bg, axis=0, return_index=True)
    sol2_bg       = sol2_bg [inds, :]
    sol1_bg       = np.vstack((sol1_bg,sol1))
    sol2_bg       = np.vstack((sol2_bg,sol2))

    # collect all the points in sol1_bg and sol2_bg, and sort them for further processing
    sol_net = np.vstack ((sol1_bg, sol2_bg))

    # start partitioning along a certain axis
    center         = np.mean (crit_points, axis=0)[:2]
    central_axis   = (crit_points[0,:2]-center)/np.linalg.norm (crit_points[0,:2]-center)

    print (f"center = {center}", flush=True)
    print (f"central axis = {central_axis}", flush=True)

    side_1  = np.empty ((0,3))
    side_2  = np.empty ((0,3))
    theta_1 = []
    theta_2 = []

    for pts in sol_net:
        direction = (pts[0:2] - center)/np.linalg.norm(pts[0:2] - center)
        clock     = np.sign(np.cross (central_axis, direction))

        if clock == 1:
            t1 = np.arccos (np.dot(direction, central_axis))
            theta_1.append (t1)
            side_1 = np.vstack ((side_1, pts))
        elif clock == -1:
            t2 = np.arccos (np.dot(direction, central_axis))
            theta_2.append (t2)
            side_2 = np.vstack ((side_2, pts))

        elif clock == 0:
            t1 = np.arccos (np.dot(direction, central_axis))
            theta_1.append (t1)
            side_1 = np.vstack ((side_1, pts))

        else:
            print ("Something's wrong.", flush=True)

    theta_1 = np.array(theta_1)
    theta_2 = np.array(theta_2)
    sorted_t1_idx = np.argsort (theta_1)
    sorted_t2_idx = np.argsort (theta_2)

    print ("Roots have been found and sorted.", flush=True)

    # WE HAVE NOW SPLIT UP THE UNSTABLE REGION INTO TWO ZONES BY USING THE CRITICAL POINT
    # NOW THAT WE HAVE TWO ZONES AS TWO HALVES:side_1 and side_2. LETS USE THEM TO FIND THE BINODAL

    # now that we have divided into two sectors, let's sort them out 
    side_1 = side_1 [sorted_t1_idx]
    side_2 = side_2 [sorted_t2_idx]

    # I now have two sides of the binodal, somewhat arbitrarily split
    # now it is time to fill in the gaps
    if len(side_1) < 100:
        side_1 = add_interpolated_rows (side_1, 5)
    if len(side_2) < 100:
        side_2 = add_interpolated_rows (side_1, 5)
    print (f"About to start refining the binodal with v2...", flush=True)
    ref_bin = refined_binodal_v2 (side_1, side_2)
    diff  = ref_bin[0] - ref_bin[1]
    norms = np.linalg.norm (diff, axis=1)
    valid = np.where (norms >= 1e-6)[0]
    ref_bin[0] = ref_bin[0][valid,:]
    ref_bin[1] = ref_bin[1][valid,:]


    # I have used fsolve to find the first set of solution points. This is going to be rather sparse. 
    # what we do now is to fill up the empty spaces using finer searches
    top_half    = np.empty ((0,3))
    bottom_half = np.empty ((0,3))
    theta_1 = []
    theta_2 = []

    for s1, s2 in zip (ref_bin[0], ref_bin[1]):
        d1 = (s1[0:2] - center)/np.linalg.norm(s1[0:2] - center)
        d2 = (s2[0:2] - center)/np.linalg.norm(s2[0:2] - center)
        clock1     = np.sign(np.cross (central_axis, d1))
        clock2     = np.sign(np.cross (central_axis, d2))

        if clock1 == 1:
            if clock2 == -1:
                pass
            else:
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
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
                print (f"This is strange. point1 = {s1}, point2 = {s2}.")
                exit ()
            t1 = np.arccos (np.dot(d1, central_axis))
            theta_2.append (t1)
            bottom_half = np.vstack ((bottom_half, s1))
            t2 = np.arccos (np.dot(d2, central_axis))
            theta_1.append (t2)
            top_half    = np.vstack ((top_half, s2))

        elif clock1 == 0:
            if clock2 == 0:
                print (f"We are at crit point.")
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

    # print (f"theta = \n{theta_1}")
    sorted_t1_idx = np.argsort (theta_1)
    top_half      = top_half [sorted_t1_idx]
    sorted_t2_idx = np.argsort (theta_2)
    bottom_half   = bottom_half[sorted_t2_idx]
    max_dist      = max_dists_on_binodal (top_half, bottom_half)

    print (f"maximum distance between two points on the binodal is {max_dist}.", flush=True)

    # this is the binodal
    ax.scatter (top_half[:,0],    top_half[:,1],    c='k', s=0.125, zorder=11)
    ax.scatter (bottom_half[:,0], bottom_half[:,1], c='k', s=0.125, zorder=11)

    if args.tl:
        for i in range (len(ref_bin[0])):
            ax.plot    ([ref_bin[0][i,0],ref_bin[1][i,0]], [ref_bin[0][i,1],ref_bin[1][i,1]], lw=0.5, ls='--', markersize=0, zorder=10)

    ff = open (args.boundary, 'w')
    ff.write ("phi_s_top|phi_p_top|phi_c_top|phi_s_bot|phi_p_bot|phi_c_bot\n")
    ff.write (f"{crit_points[0,0]}|{crit_points[0,1]}|{1-crit_points[0,0]-crit_points[0,1]}|{crit_points[1,0]}|{crit_points[1,1]}|{1-crit_points[1,0]-crit_points[1,1]}\n")
    for i in range (len(ref_bin[0])):
        ff.write (f"{ref_bin[0][i][0]}|{ref_bin[0][i][1]}|{ref_bin[0][i][2]}|{ref_bin[1][i][0]}|{ref_bin[1][i][1]}|{ref_bin[1][i][2]}\n")

    return




if __name__=="__main__":

    start = time.time()

    ###########################
    N = args.N
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc
    chi_ac = args.chi_ac
    nproc  = args.nproc
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

    binodal_plotter_prototype (fig, ax, dumpfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc, crits)
    print ("Done with binodal plotting!", flush=True)

    ax.plot (crits[:,0], crits[:,1], c='darkred', zorder=12, lw=0.5)

    # Set axis limits
    ax.set_xlim (0, 1)
    ax.set_ylim (0, 1)

    ax.grid()

    fig.savefig  (args.img, dpi=1200)
    print ("Completed heat map computation.", flush=True)
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.", flush=True)
