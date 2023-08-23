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

################################

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
            root = fsolve (mu_equations, [bincloser[gidx, 1], binfurther[gidx, 0], binfurther[gidx,1]])
            p1   = np.array ([phi_a, root[0], 1-phi_a-root[0]])
            p2   = np.array ([root[1], root[2], 1-root[1]-root[2]])
            if (np.abs ( np.array (mu_equations (root))) > 1e-6).any ():
                continue

            elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                continue

            elif (np.linalg.norm(sol1_bg - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2_bg - p2, axis=-1)<1e-6).any():
                continue

            elif np.linalg.norm (p1-p2) < 1e-6:
                continue

            else:
                # print (f"HIT!")
                sol1_bg = np.vstack ((sol1_bg, p1))
                sol2_bg = np.vstack ((sol2_bg, p2))
                break

    return [sol1_bg, sol2_bg]


def binodal_plotter (fig, ax, dumpfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc, crit_points):

    try:
        df = pd.read_csv (dumpfile, sep='\s+', engine="python", skiprows=1, names=["dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2"])
        df = df.loc[df["dmu"]<5]
        print (df)
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
            eq1 = mu_a(phi_a[idx], phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi_a[idx], phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi_a[idx], phi[0]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_b[idx], phi_an[idx], phi_bn[idx]])
        binodal_closer [idx, :] = np.array([phi_a[idx], root[0], 1-root[0]-phi_a[idx]])
        binodal_further[idx, :] = np.array([root[1], root[2], 1-root[1]-root[2]])
        # print (f"phi_b[{idx}] = {phi_b[idx]}")

        # if the roots are "bad" roots, just write them out as bad
        if ( np.abs(np.array(mu_equations(root))) > 1e-6 ).any():
            # print ([phi_a[idx], root[0], 1-phi_a[idx]-root[0]], [root[1], root[2], 1-root[1]-root[2]],", error = ",mu_equations(root))
            bad_idx.append (idx)
        else:
            fa = [phi_a[idx], root[1]]
            fb = [root[0], root[2]]
            fc = [1-root[0]-phi_a[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_a[idx], 1-root[0]-phi_a[idx]])
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
    # sol1_bg, inds = np.unique (sol1_bg, axis=0, return_index=True)
    # sol2_bg       = sol2_bg [inds, :]
    sol1_bg       = np.vstack((sol1_bg,sol1))
    sol2_bg       = np.vstack((sol2_bg,sol2))

    ref_bin = [sol1_bg, sol2_bg]

    # print (f"sol1_bg = {sol1_bg}")
    # print (f"sol2_bg = {sol2_bg}")
    # this is the binodal
    # ax.scatter (ref_bin[0][:,0], ref_bin[0][:,1], s=0.125, zorder=11, c="darkred")
    # ax.scatter (ref_bin[1][:,0], ref_bin[1][:,1], s=0.125, zorder=11, c="darkred" )

    if args.tl:
        for i in range (len(ref_bin[0])):
            line = np.linspace (ref_bin[0][i,0:2], ref_bin[1][i,0:2], 100)
            if (stab_crit (line[:,0], line[:,1], chi_ab, chi_bc, chi_ac) < 0).any():
                ax.plot    ([ref_bin[0][i,0],ref_bin[1][i,0]], [ref_bin[0][i,1],ref_bin[1][i,1]], lw=0.5, ls='--', markersize=0, zorder=10)
                ax.scatter (ref_bin[0][i,0], ref_bin[0][i,1], s=0.125, zorder=11, c="darkred")
                ax.scatter (ref_bin[1][i,0], ref_bin[1][i,1], s=0.125, zorder=11, c="darkred" )
            else:
                pass

    sorter = ref_bin[0][:,1].argsort()
    ref_bin[0] = ref_bin[0][sorter]
    ref_bin[1] = ref_bin[1][sorter]

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
    nproc  = args.nproc
    dumpfile = args.dumpfile
    ############################

    lsize = 3
    font  = {'color':  'black', 'weight': 'normal', 'size': lsize}
    fig   = plt.figure(num=1, figsize=(3,3))
    ax    = plt.axes ()

    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4* N * (1 - 2* phi_s * chi_sc + 2 * phi_s ** 2 * chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    (-1 + 2 * phi_s * chi_sc + N * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc:  1 / (2*N * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: 1 - 2 * phi_s * chi_sc + N * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )


    roots_up, roots_down = find_crit_point (N, chi_ac, chi_ab, chi_bc)
    ax.scatter (roots_up[:,0]  , roots_up[:,1],   color='k', edgecolors='steelblue', s=0.1, zorder=11)
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
    p_a_space = np.linspace (0.001, 1-0.001, 1000)
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

    binodal_plotter (fig, ax, dumpfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc, crits)
    print ("Done with binodal plotting! Saving image...", flush=True)

    # Set axis limits
    ax.set_xlim (0, 1)
    ax.set_ylim (0, 1)

    ax.grid ()

    fig.savefig  (args.img, dpi=1200)
    print (f"Completed heat map computation.", flush=True)
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.", flush=True)

