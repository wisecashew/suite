import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import argparse
import time
import matplotlib.pyplot as plt
import mpltern

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
    df = pd.read_csv ("comparing_mu_hold.txt", sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
        "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])

    df = df.loc[df["dmu"]<0.01]

    # print (df)

    va = 1
    vb = 32
    N  = 32
    vc = 1

    chi_ab = -1
    chi_ac = 0
    chi_bc = -10

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )


    phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
    phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
    phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

    phi_a1 = np.zeros (phi_a.shape); phi_a2 = np.zeros (phi_an.shape);
    phi_b1 = np.zeros (phi_b.shape); phi_b2 = np.zeros (phi_bn.shape);
    phi_c1 = np.zeros (phi_c.shape); phi_c2 = np.zeros (phi_cn.shape);

    # now it is time to find the best guesses. Check each value of phi_a, phi_an and see if they lie inside or outside the spinodal. 
    # if both lie outside, great. They are good guesses. If one of them doesn't extend the line in the direction where it 
    # doesn't until you hit the boundary of the spinodal. If both of them lie inside, kill both.

    # stability criterion
    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    for i in range (len (phi_a)):
        p1 = np.array ([phi_a[i] , phi_b[i] , phi_c[i]] )
        p2 = np.array ([phi_an[i], phi_bn[i], phi_cn[i]])

        if np.sign (stab_crit(p1[0], p1[1], chi_ab, chi_bc, chi_ac)) == 1:
            r = p2 - p1
            extension = lambda L: stab_crit ((p1 + L*r)[0],(p1 + L*r)[1], chi_ab, chi_bc, chi_ac)
            scale = 2
            root      = fsolve (extension, scale)
            p2        = p1 + root * r
            while (np.sum(p2)>1):
                scale = scale/1.1
                root      = fsolve (extension, scale)
                p2        = p1 + root * r
            
        
        else:
            r = p1 - p2
            extension = lambda L: stab_crit ((p1 + L*r)[0],(p1 + L*r)[1], chi_ab, chi_bc, chi_ac)
            scale     = 2
            root      = fsolve (extension, scale)
            p1        = p2 + root * r
            while (np.sum(p2)>1):
                scale = scale/1.1
                root      = fsolve (extension, scale)
                p2        = p1 + root * r

        phi_a1[i] = p1[0]; phi_a2[i] = p2[0]; 
        phi_b1[i] = p1[1]; phi_b2[i] = p2[1];
        phi_c1[i] = p1[2]; phi_c2[i] = p2[2];

    
    for idx in range (len(phi_a)):

        def mu_equations (phi):
            eq1 = mu_a(phi[0], phi_b1[idx]) - mu_a(phi[1], phi[2])
            eq2 = mu_b(phi[0], phi_b1[idx]) - mu_b(phi[1], phi[2])
            eq3 = mu_c(phi[0], phi_b1[idx]) - mu_c(phi[1], phi[2])

            return [eq1, eq2, eq3]

        root = fsolve (mu_equations, [phi_a1[idx], phi_a2[idx], phi_b2[idx]])

        if ( np.array(mu_equations(root)) > 1e-6).any():
            pass
        else:
            fa = [root[0], root[1]]
            fb = [phi_b1[idx], root[2]]
            fc = [1-root[0]-phi_b1[idx], 1-root[1]-root[2]]
            p1 = np.array([root[0], phi_b1[idx], 1-root[0]-phi_b1[idx]])
            p2 = np.array([root[1], root[2], 1-root[1]-root[2]])

            if np.linalg.norm (p1-p2) > 0.05:
                ax.scatter (fa, fb, fc, s=1, c='steelblue')
                ax.plot    (fa, fb, fc, lw=0.5, ls='--', markersize=0)

    ax.grid ()
    positions = ['tick1', 'tick2']
    for position in positions:
        ax.taxis.set_ticks_position(position)
        ax.laxis.set_ticks_position(position)
        ax.raxis.set_ticks_position(position)

    ax.set_tlabel('Vol. frac. A')
    ax.set_llabel('Vol. frac. B')
    ax.set_rlabel('Vol. frac. C')

    fig.savefig ("binodal_with_lines", dpi=1200)
    stop = time.time()
    print (f"Elapsed time is {stop-start} seconds.")