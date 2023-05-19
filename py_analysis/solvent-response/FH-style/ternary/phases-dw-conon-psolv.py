#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

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
import mpltern
import sys
np.set_printoptions(threshold=sys.maxsize)


if __name__=="__main__":


    lsize = 3
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(figsize=(5,5))
    ax  = fig.add_subplot (projection="ternary")
    
    N = 100

    discriminant = lambda phi_a, chi_ab, chi_bc, chi_ac: -4* N * (1 - 2* phi_a * chi_ac + 2 * phi_a ** 2 * chi_ac) * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) **2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) + \
    (-1 + 2 * phi_a * chi_ac + N * (1 - 2*chi_bc - phi_a * (chi_ab ** 2 + chi_ac **2 - 2*chi_ac*chi_bc + (chi_bc -2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) + phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) ) ** 2

    denom  = lambda phi_a, chi_ab, chi_bc, chi_ac:  1 / (2*N * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2*chi_ab * (chi_ac + chi_bc) ) ) )

    prefac = lambda phi_a, chi_ab, chi_bc, chi_ac: 1 - 2 * phi_a * chi_ac + N * ( -1 + 2 * chi_bc + phi_a * (chi_ab ** 2 + chi_ac ** 2 - 2*chi_ac * chi_bc + (chi_bc - 2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) - phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) )

    roots  = lambda phi_a, chi_ab, chi_bc, chi_ac: np.sort(np.array([denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) + np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) ) , denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) - np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) ) ]))


    meshsize            = 1000
    phi_a               = np.linspace (0.01, 1-0.01, meshsize)
    chi_ab              = -1
    chi_bc              = -1
    chi_ac              = -12

    tern = roots (phi_a, chi_ab, chi_bc, chi_ac)
    r1   = tern[0]
    r2   = tern[1]

    to_keep = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0) * (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)

    r1    = r1[to_keep]
    r2    = r2[to_keep]
    phi_a = phi_a[to_keep]

    # print (r1)
    
    # Plot the points
    ax.scatter(phi_a, r1, 1-phi_a-r1, color='steelblue', s=1)
    ax.scatter(phi_a, r2, 1-phi_a-r2, color='coral', s=1)

    ax.set_tlabel('Vol. frac. A')
    ax.set_llabel('Vol. frac. B')
    ax.set_rlabel('Vol. frac. C')

    # Set axis limits
    ax.set_tlim(0, 1)
    ax.set_llim(0, 1)
    ax.set_rlim(0, 1)
    # ax.ticks(axis='lbr', multiple=5, linewidth=1, offset=0.025)

    positions = ['tick1', 'tick2']
    for position in positions:
        ax.taxis.set_ticks_position(position)
        ax.laxis.set_ticks_position(position)
        ax.raxis.set_ticks_position(position)

    # Add gridlines
    # gridlines_per_axis = 100
    # ax.grid(linewidth=0.5, color='gray', linestyle='dotted', multiple=gridlines_per_axis)
    ax.grid()

    plt.savefig ("ternary", dpi=1200)
    """
    print ("Computation completed.")

    # print ("solvent_phi_range = ",solvent_phi_range)
    """
