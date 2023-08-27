#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
from scipy.spatial import ConvexHull
import mpltern
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
parser.add_argument('--skelfile',     dest='skelfile',  type=str,            action='store',  help="name of file where the skeleton was dumped.")
parser.add_argument('--bin-boundary', dest='boundary',  type=str,            action='store',  help="name of file where you will dump out your solution for the binodal")
parser.add_argument('--tielines',     dest='tl',        action='store_true', default=False,   help="Option to include if you want to see all tie-lines.")
parser.add_argument('--image',        dest='img',       type=str,            action='store',  help="name of image generated.")
args = parser.parse_args()


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
            elif (np.linalg.norm(sol1_bg - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2_bg - p2, axis=-1)<1e-6).any():
                continue
            
            elif np.linalg.norm (p1-p2) < 1e-6:
                continue

            else:
                sol1_bg = np.vstack ((sol1_bg, p1))
                sol2_bg = np.vstack ((sol2_bg, p2))
                break

    return (sol1_bg, sol2_bg)


# binodal_plotter, as the name suggests, plots binodals
def binodal_plotter (fig, ax, fig2, ax2, skelfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc):

    df = pd.read_csv (skelfile, sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])
    print (df["dmu"])
    df = df.loc[df["dmu"]<1]

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
    f = open (args.boundary, 'w')
    print ("Start processing the dumpfile and find roots.", flush=True)
    for idx in range (len(phi_a)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]])
        binodal_closer [idx, :] = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
        binodal_further[idx, :] = np.array([root[1], root[2], 1-root[1]-root[2]])

        # if the roots are "bad" roots, just write them out as bad
        if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
            bad_idx.append (idx)
            f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")
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
                    f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")

                elif stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) < 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) < 0:
                     bad_idx.append (idx)
                     f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")

                elif np.isnan(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isnan(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)
                     f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")

                elif np.isinf(stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac)) or np.isinf(stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac)):
                     bad_idx.append (idx)
                     f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")


                else:
                    good_idx.append (idx)
                    sol1 = np.vstack ((sol1,p1))
                    sol2 = np.vstack ((sol2,p2))
                    f.write (f"idx = {idx}: good  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")

    f.close ()
    print ("Everything has been written out. Start doing a bigger processing - with parallelization.", flush=True)

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
    
    # NEXT, only keep things on the surface of the binodal, none of the inside stuff
    # find the convex hull in order to do the above

    points1 = np.vstack ((sol1[:, 0:2], sol1_bg[:, 0:2]))
    points2 = np.vstack ((sol2[:, 0:2], sol2_bg[:, 0:2]))
    
    if args.tl:
        for i in range (len(sol1_bg)):
            ax.plot    ([sol1_bg[i][0], sol2_bg[i][0]], [sol1_bg[i][1],sol2_bg[i][1]], [1-sol1_bg[i][0]-sol1_bg[i][1], 1-sol2_bg[i][0]-sol2_bg[i][1]], lw=0.1, ls='--', markersize=0, c='lightgray')

    f = open (args.boundary, 'w')

    for i in range(sol1.shape[0]):
        f.write (f"{sol1[i,0]}, {sol1[i,1]}, {sol1[i,2]}, {sol2[i,0]}, {sol2[i,1]}, {sol2[i,2]}\n")
    for i in range(sol1_bg.shape[0]):
        f.write (f"{sol1_bg[i,0]}, {sol1_bg[i,1]}, {sol1_bg[i,2]}, {sol2_bg[i,0]}, {sol2_bg[i,1]}, {sol2_bg[i,2]}\n")

    f.close ()


    return



if __name__=="__main__":

    start = time.time()

    ###########################
    N = args.N
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc
    chi_ac = args.chi_ac
    nproc  = args.nproc
    skelfile = args.skelfile
    ############################

    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure (num=1, figsize=(6, 6))
    ax  = fig.add_subplot ()

    fig2  = plt.figure (num=2)
    ax2   = plt.axes   ()


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
    p_a  = p_a[to_keep]
    p_b  = p_b[to_keep]

    vmax = np.max (vals)
    vmin = np.min (vals)

    if np.sign (vmax) == np.sign (vmin):
        if np.sign (vmax) >=0:
            vmin = -vmax
        else:
            vmax = -vmin
            
    print (f"vmin = {vmin}, vmax = {vmax}")
    
    norm = colors.SymLogNorm (0.001, vmin=vmin, vmax=vmax)
    cols = cm.bwr (norm (vals))
    
    # Plot the points
    p_c = 1 - p_a - p_b
    ax.scatter  (p_a, p_b, s=0.01, color=cols)
    print ("Painted the diagram!", flush=True)
    
    ax.set_xlabel ("$\\phi _{A}$")
    ax.set_ylabel ("$\\phi _{B}$")
    
    print ("We have plotted the spinodal region!\n\n")
    
    print ("###########################################################\n\n")
    
    print ("Start binodal plotting...")
    
    va = 1
    vb = N
    vc = 1
    
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )
    
    binodal_plotter (fig, ax, fig2, ax2, skelfile, nproc, chi_ab, chi_bc, chi_ac, va, vb, vc)
    print   ("Done with binodal plotting!")
    ax.grid ()
    
    fig.savefig  (args.img, dpi=1200)
    
    print ("Completed heat map computation.")
    stop = time.time ()
    print (f"Time for computation is {stop-start} seconds.")
    
    