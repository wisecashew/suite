import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import argparse
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import mpltern
import argparse
from matplotlib.collections import LineCollection

parser = argparse.ArgumentParser(description="Create a ternary spinodal diagram.")
parser.add_argument('--dumpfile', dest='dumpfile', type=str, action='store', help="name of dumpfile.")
parser.add_argument('--image', dest='img', type=str, action='store', help="name of image generated.")
args = parser.parse_args() 


if __name__=="__main__":

    start = time.time()
    lsize = 3
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(5,5))
    ax  = fig.add_subplot (projection="ternary")

    # binodal plotter
    df = pd.read_csv (args.dumpfile, sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])

    df = df.loc[df["dmu"]<1]

    va = 1
    vb = 32
    vc = 1

    chi_ab = -1
    chi_ac = 0
    chi_bc = -10
    N      = 32

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

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    sol1 = np.empty((0,3))
    sol2 = np.empty((0,3))
    bad_idx  = []
    good_idx = []
    f = open ("root-finding.out", 'w')
    for idx in range (len(phi_a)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]])
        binodal_closer [idx, :] = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
        binodal_further[idx, :] = np.array([root[1], root[2], 1-root[1]-root[2]])
        # print (f"phi_b[{idx}] = {phi_b[idx]}")

        if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
            bad_idx.append (idx)
            f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")
        else:
            fa = [root[0], root[1]]
            fb = [phi_b[idx], root[2]]
            fc = [1-root[0]-phi_b[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if np.linalg.norm (p1-p2) > 0.05:
                if (np.linalg.norm(sol1 - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2 - p2, axis=-1)<1e-6).any():
                    bad_idx.append (idx)
                    f.write (f"idx = {idx}: bad  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")
                    pass
                else:
                    good_idx.append (idx)
                    sol1 = np.vstack ((sol1,p1))
                    sol2 = np.vstack ((sol2,p2))
                    f.write (f"idx = {idx}: good  points = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}, {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}, roots: {root[0], phi_b[idx], 1-root[0]-phi_b[idx]}, {root[1], root[2], 1-root[1]-root[2]}\n")

    # now, time to slowly convert the "bad" points to "good" points
    sol1_bg = np.empty((0,3))
    sol2_bg = np.empty((0,3))
    print (f"Size of bad_idx = {len(bad_idx)}")
    for i, idx in enumerate (bad_idx):
        print (f"at bad_idx = {idx}...", flush=True)
        phi_a = binodal_closer[idx,0]
        def mu_equations (phi):
            eq1 = mu_a(phi_a, phi[0]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi_a, phi[0]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi_a, phi[0]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        # now go through every single good index 
        for g, gidx in enumerate (good_idx):
            root = fsolve (mu_equations, [init[gidx,1], binodal_closer[gidx,0], binodal_closer[gidx,1]])
            p1 = np.array ([phi_a, root[0], 1-phi_a-root[0]])
            p2 = np.array ([root[1], root[2], 1-root[1]-root[2]])
            if (np.abs (np.array (mu_equations (root))) > 1e-6).any ():
                continue

            elif (np.linalg.norm(sol1_bg - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2_bg - p2, axis=-1)<1e-6).any():
                continue

            elif np.linalg.norm (p1-p2) < 1e-6: # or ( stab_crit (p1[0], p1[1], chi_ab, chi_bc, chi_ac) > 0 or stab_crit (p2[0], p2[1], chi_ab, chi_bc, chi_ac) > 0 )):
                continue

            else:
                sol1_bg = np.vstack((sol1_bg, p1))
                sol2_bg = np.vstack((sol2_bg, p2))
                # print (f"p1 = {p1}")
                # print (f"p2 = {p2}")
                break



    ax.scatter (sol1_bg[:,0], sol1_bg[:,1], sol1_bg[:,2], s=1, c='steelblue')
    ax.scatter (sol2_bg[:,0], sol2_bg[:,1], sol2_bg[:,2], s=1, c='steelblue')
    ax.scatter (sol1[:,0], sol1[:,1], sol1[:,2], s=1, c='steelblue')
    ax.scatter (sol2[:,0], sol2[:,1], sol2[:,2], s=1, c='steelblue')
    solnet = np.vstack ((sol1, sol2))

    for i in range (sol1.shape[0]):
        ax.plot    ([sol1[i,0],sol2[i,0]], [sol1[i,1],sol2[i,1]], [sol1[i,2],sol2[i,2]], lw=0.1, ls='--', markersize=0)
    for i in range (sol1_bg.shape[0]):
        ax.plot    ([sol1_bg[i,0],sol2_bg[i,0]], [sol1_bg[i,1],sol2_bg[i,1]], [sol1_bg[i,2],sol2_bg[i,2]], lw=0.1, ls='--', markersize=0)

    ax.grid ()
    positions = ['tick1', 'tick2']
    for position in positions:
        ax.taxis.set_ticks_position(position)
        ax.laxis.set_ticks_position(position)
        ax.raxis.set_ticks_position(position)

    ax.set_tlabel('Vol. frac. A')
    ax.set_llabel('Vol. frac. B')
    ax.set_rlabel('Vol. frac. C')

    fig.savefig (args.img, dpi=1200)
    stop = time.time()

    f.close ()

    print (f"Elapsed time is {stop-start} seconds.", flush=True)
