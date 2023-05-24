#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.patheffects as pe 
from scipy.optimize import fsolve
from scipy.optimize import root
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator


if __name__=="__main__":


    lsize = 8
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(figsize=(2.5,2.0), constrained_layout=True)
    fig.tight_layout()
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=lsize, pad=5)
    ax.tick_params(axis='y', labelsize=lsize)

    
    chi = lambda E, T: E/T
    N   = 100

    frac_a = lambda r, phi_b: (1-phi_b)*r/(r+1)
    D      = lambda phi_a, phi_b, E_ab, E_bc, E_ac, T: (1/(N*phi_b) + 1/(1-phi_a-phi_b) - 2*chi(E_bc, T)) * (1/(phi_a) + 1/(1-phi_a-phi_b) - 2*chi(E_ab, T)) - (1/(1-phi_a-phi_b) + chi(E_ab,T) - chi(E_bc,T) - chi(E_ac,T)) ** 2

    # I will now attempt to create Fig. 1 from Dudowicz's paper
    # this is for cosolvency 
    E_ab  = 300
    E_bc  = 300
    ratio = 1
    T_lower = []
    T_upper = []
    my_phi  = []

    phi_range = np.arange (0.4, 0.5, 0.02)
    seeds     = [10000]


    for E_ac in [0, 100, 200, 300, 400]:

        T_lower.clear ()
        T_upper.clear ()
        my_phi.clear  ()

        for p_b in phi_range:
            
            D_specified = lambda T: D (frac_a (ratio, p_b), p_b, E_ab, E_bc, E_ac, T)
            # print (f"phi_a = {frac_a (ratio, p_b)}")
            roots       = np.array ([])
            roots_error = np.array ([])

            for s in seeds:
                roots = np.hstack ( ( roots, root (D_specified, 1000, method='hybr', tol=1e-9) ) )
                roots_error = np.hstack ( ( roots_error, abs ( D_specified(roots[-1]) ) ) )

            hold        = roots > 0
            roots       = roots[hold]
            roots_error = roots_error [hold]

            hold        = roots_error < 1e-9
            roots       = roots[hold]
            roots_error = roots_error [hold]
            
            try:
                T_lower.append (np.min (roots))
                T_upper.append (np.max (roots))
                my_phi.append  (p_b)
            except ValueError:
                print ("Root finding was unstable.")

        ax.plot (my_phi[::1], T_lower[::1], ls='-', lw=1, zorder=10, solid_capstyle='round') # , path_effects=[pe.Stroke(linewidth=3.5, foreground='k'), pe.Normal()]) # , clip_on=False)
        ax.plot (my_phi[::1], T_upper[::1], ls='-', lw=1, zorder=10, solid_capstyle='round', label=f"$\\chi _{{ac}} = {E_ac}/T$") # , path_effects=[pe.Stroke(linewidth=3.5, foreground='k'), pe.Normal()]) #, clip_on=False)



    ax.minorticks_on ()
    ax.set_xticks (np.linspace (0, 0.5, 6))
    ax.legend(fontsize="xx-small")
    # ax.set_yscale ("log")
    # ax.set_yticks ([0.1, 1, 10, 50 ])
    # ax.set_yticklabels (ax.get_yticks(), fontdict=font) # , weight='bold')
    # if pv == 1.0:
    #     ax.set_xticklabels (ax.get_xticks(), fontdict=font, ) # , weight='bold')
    #     plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    # else:
    #     ax.set_xticklabels([])
    
    # ax.set_ylim (0.1, 50)
    ax.set_xlim (0, 0.5)
    # ax.legend(loc="upper right", fontsize=4, frameon=False, ncol=2)
    plt.savefig ("plots.png", dpi=1200) # bbox_inches='tight', dpi=1200)
