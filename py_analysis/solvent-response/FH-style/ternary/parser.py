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
    ax.taxis.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.laxis.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.raxis.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    # binodal plotter
    df = pd.read_csv (args.dumpfile, sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])

    df = df.loc[df["dmu"]<0.1]
    df = df.loc[df["phi_b1"]<0.01]
    # df = df.loc[df["phi_b1"]<0.42]
    df = df.loc[df["phi_b2"]>0.4]
    # df.to_csv ("selective-output.csv")
    # print (df[["dmu", "phi_a1", "phi_b1", "phi_a2", "phi_b2"]])
    # exit ()

    va = 1
    vb = 32
    vc = 1

    chi_ab = -1
    chi_ac = 0
    chi_bc = -10
    N      = 32

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2


    phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
    phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
    phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    sol1 = np.empty((0,3))
    sol2 = np.empty((0,3))

    for idx in range (len(phi_a))[0:10]:

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a[idx], phi_an[idx], phi_bn[idx]])
        p1 = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
        p2 = np.array([root[1], root[2], 1-root[1]-root[2]])
        print (f"guess1 = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}")
        print (f"guess2 = {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}")
        print (f"root1 = {p1}")
        print (f"root2 = {p2}")
        print ("##########################\n")
        if ( np.abs(np.array(mu_equations(root))) > 1e-6).any():
            pass
        else:
            fa = [root[0], root[1]]
            fb = [phi_b[idx], root[2]]
            fc = [1-root[0]-phi_b[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b[idx], 1-root[0]-phi_b[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if np.linalg.norm (p1-p2) > 0.05:
                if (np.linalg.norm(sol1 - p1, axis=-1)<1e-6).any() or (np.linalg.norm(sol2 - p2, axis=-1)<1e-6).any():
                    pass
                else:
                    sol1 = np.vstack ((sol1,p1))
                    sol2 = np.vstack ((sol2,p2))
                    print (f"Good roots found!")
                    print (f"g1 = {phi_a[idx], phi_b[idx], 1-phi_a[idx]-phi_b[idx]}")
                    print (f"g2 = {phi_an[idx], phi_bn[idx], 1-phi_an[idx]-phi_bn[idx]}")
                    print (f"p1 = {p1}")
                    print (f"p2 = {p2}")
                    print (f"roots = {root}")
                    print (f"delta mu = {mu_equations(root)}")
    exit ()
    ax.scatter (sol1[:,0], sol1[:,1], sol1[:,2], s=1, c='steelblue')
    ax.scatter (sol2[:,0], sol2[:,1], sol2[:,2], s=1, c='steelblue')
    solnet = np.vstack ((sol1, sol2))

    for i in range (sol1.shape[0]):
        ax.plot    ([sol1[i,0],sol2[i,0]], [sol1[i,1],sol2[i,1]], [sol1[i,2],sol2[i,2]], lw=0.1, ls='--', markersize=0)

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
    print (f"Elapsed time is {stop-start} seconds.")
