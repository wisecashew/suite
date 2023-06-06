import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import mpltern
import sys
import argparse
np.set_printoptions(threshold=sys.maxsize)

import argparse 
parser = argparse.ArgumentParser(description="Create two spinodal diagram: one with only shows the edges (ternary) and another one which paints the ternary plot.")
parser.add_argument('--chiac', metavar='chi_ac', dest='chi_ac', type=float, action='store', help='enter A-C exchange parameter.')
parser.add_argument('--chiab', metavar='chi_ab', dest='chi_ab', type=float, action='store', help='enter A-B exchange parameter.')
parser.add_argument('--chibc', metavar='chi_bc', dest='chi_bc', type=float, action='store', help='enter B-C exchange parameter.')
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
args = parser.parse_args() 


if __name__=="__main__":


    lsize = 3
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(5,5))
    ax  = fig.add_subplot (projection="ternary")
    
    N = args.N

    # FIND PHI_B GIVEN PHI_A

    discriminant = lambda phi_a, chi_ab, chi_bc, chi_ac: -4* N * (1 - 2* phi_a * chi_ac + 2 * phi_a ** 2 * chi_ac) * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) **2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) + \
    (-1 + 2 * phi_a * chi_ac + N * (1 - 2*chi_bc - phi_a * (chi_ab ** 2 + chi_ac **2 - 2*chi_ac*chi_bc + (chi_bc -2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) + phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) ) ** 2

    denom  = lambda phi_a, chi_ab, chi_bc, chi_ac:  1 / (2*N * (2*chi_bc + phi_a * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2*chi_ab * (chi_ac + chi_bc) ) ) )

    prefac = lambda phi_a, chi_ab, chi_bc, chi_ac: 1 - 2 * phi_a * chi_ac + N * ( -1 + 2 * chi_bc + phi_a * (chi_ab ** 2 + chi_ac ** 2 - 2*chi_ac * chi_bc + (chi_bc - 2) * chi_bc - 2 * chi_ab * (-1 + chi_ac + chi_bc) ) - phi_a ** 2 * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) )

    root_up  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) + np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )
    root_lo  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) - np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )

    meshsize            = 1000
    phi_a               = np.linspace (0.001, 1-0.001, meshsize*10)
    chi_ab              = args.chi_ab
    chi_bc              = args.chi_bc
    chi_ac              = args.chi_ac

    r1 = root_up (phi_a, chi_ab, chi_bc, chi_ac)
    r2 = root_lo (phi_a, chi_ab, chi_bc, chi_ac)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 2) * (r2 >= 0)
    r2 = r2[to_keep_2]

    
    # Plot the points
    ax.scatter(phi_a[to_keep_1], r1, 1-phi_a[to_keep_1]-r1, color='maroon', s=1)
    ax.scatter(phi_a[to_keep_2], r2, 1-phi_a[to_keep_2]-r2, color='salmon', s=1)

    # FIND PHI_A GIVEN PHI_B

    discriminant = lambda phi_b, chi_ab, chi_bc, chi_ac: -4 * (1 + 2 * N * phi_b ** 2 * chi_bc + phi_b * (-1 + N - 2 * N * chi_bc) ) * ( 2 * chi_ac + N * phi_b * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) ) +\
    (2 * (-1 + phi_b) * chi_ac + N * phi_b * ( (-1 + phi_b) * chi_ab ** 2 + (-1 + phi_b) * chi_ac ** 2 - 2 * (-1 + phi_b) * chi_ac * chi_bc + chi_bc * (2 + (-1 + phi_b) * chi_bc ) - 2 * chi_ab * (1 + (-1 + phi_b) * chi_ac + (-1 + phi_b) * chi_bc ) ) ) ** 2 

    denom    = lambda phi_b, chi_ab, chi_bc, chi_ac:  1 / (4 * chi_ac + 2 * N * phi_b * (chi_ab ** 2 + (chi_ac - chi_bc) ** 2 - 2 * chi_ab * (chi_ac + chi_bc) ) )

    prefac   = lambda phi_b, chi_ab, chi_bc, chi_ac: -2 * (-1 + phi_b) * chi_ac + N * phi_b * ( - ( (-1 + phi_b) * chi_ab ** 2) - (-1 + phi_b) * chi_ac ** 2 + 2 * (-1 + phi_b) * chi_ac * chi_bc + 2 * chi_ab * (1 + (-1 + phi_b) * chi_ac + (-1 + phi_b) * chi_bc ) + chi_bc * (-2 + chi_bc - phi_b * chi_bc) )

    root_up  = lambda phi_b, chi_ab, chi_bc, chi_ac: denom (phi_b, chi_ab, chi_bc, chi_ac) * ( prefac (phi_b, chi_ab, chi_bc, chi_ac) + np.sqrt(discriminant (phi_b, chi_ab, chi_bc, chi_ac) ) )
    root_lo  = lambda phi_b, chi_ab, chi_bc, chi_ac: denom (phi_b, chi_ab, chi_bc, chi_ac) * ( prefac (phi_b, chi_ab, chi_bc, chi_ac) - np.sqrt(discriminant (phi_b, chi_ab, chi_bc, chi_ac) ) )

    phi_b    = np.linspace (0.001, 1-0.001, meshsize*10)

    r1 = root_up (phi_b, chi_ab, chi_bc, chi_ac)
    r2 = root_lo (phi_b, chi_ab, chi_bc, chi_ac)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 2) * (r2 >= 0)
    r2 = r2[to_keep_2]

    
    # Plot the points
    ax.scatter (r1, phi_b[to_keep_1], 1-phi_b[to_keep_1]-r1, color='darkorange', s=1)
    ax.scatter (r2, phi_b[to_keep_2], 1-phi_b[to_keep_2]-r2, color='moccasin',   s=1)
    # ax.scatter(phi_b[to_keep_1], r1, 1-phi_b[to_keep_1]-r1, color='darkorange', s=1)
    # ax.scatter(phi_b[to_keep_2], r2, 1-phi_b[to_keep_2]-r2, color='moccasin',   s=1)


    discriminant = lambda phi_a, chi_ab, chi_bc, chi_ac: -4 * N * (phi_a + N * (-1 + phi_a) * (-1 + 2 * chi_ab * phi_a) ) * ( (chi_ab - chi_ac) ** 2 * phi_a + chi_bc ** 2 * phi_a - 2 * chi_bc * (-1 + chi_ab * phi_a + chi_ac * phi_a) ) +\
    (-1 + 2 * chi_ac * phi_a + N * (1 + chi_ac ** 2 * phi_a + 2 * chi_ab * (-1 + chi_ac * (-1 + phi_a) ) * phi_a - chi_ab ** 2 * (phi_a - 1) * phi_a - chi_bc ** 2 * (-1 + phi_a) * phi_a - chi_ac ** 2 * phi_a ** 2 + 2 * chi_bc * (-1 + phi_a) * (-1 + chi_ab * phi_a + chi_ac * phi_a) ) ) ** 2

    denom        = lambda phi_a, chi_ab, chi_bc, chi_ac:  1 / (-2 * N * ( (chi_ab - chi_ac) ** 2 * phi_a + chi_bc ** 2 * phi_a - 2 * chi_bc * (-1 + chi_ab * phi_a + chi_ac * phi_a) ) )

    prefac       = lambda phi_a, chi_ab, chi_bc, chi_ac: 1 - 2 * chi_ac * phi_a + N * (-1 - chi_ac ** 2 * phi_a + chi_ab ** 2 * (-1 + phi_a) * phi_a + chi_bc ** 2 * (-1 + phi_a) * phi_a + chi_ac ** 2 * phi_a ** 2 + 2 * chi_ab * phi_a * (1 + chi_ac - chi_ac * phi_a) - 2 * chi_bc * (-1 + phi_a) * (-1 + chi_ab * phi_a + chi_ac * phi_a) )

    root_up  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) + np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )
    root_lo  = lambda phi_a, chi_ab, chi_bc, chi_ac: denom (phi_a, chi_ab, chi_bc, chi_ac) * ( prefac (phi_a, chi_ab, chi_bc, chi_ac) - np.sqrt(discriminant (phi_a, chi_ab, chi_bc, chi_ac) ) )

    phi_a    = np.linspace (0.001, 1-0.001, meshsize*10)

    r1 = root_up (phi_a, chi_ab, chi_bc, chi_ac)
    r2 = root_lo (phi_a, chi_ab, chi_bc, chi_ac)

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 2) * (r2 >= 0)
    r2 = r2[to_keep_2]

    
    # Plot the points
    ax.scatter(phi_a[to_keep_1], 1-phi_b[to_keep_1]-r1, r1, color='forestgreen', s=1)
    ax.scatter(phi_a[to_keep_2], 1-phi_b[to_keep_2]-r2, r2, color='springgreen', s=1)


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
    # ax.laxis.set_major_locator (MultipleLocator (6))
    ax.grid()

    plt.savefig (f"edges_{chi_ab}_{chi_bc}_{chi_ac}.png", dpi=1200)
    
    print ("Completed spinodal boundary plotting!\n\n")
    print (" ########################################### \n\n")
    print ("Start painting the spinodal region...")

    
    lsize = 3
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=2, figsize=(5,5))
    ax  = fig.add_subplot (projection="ternary")

    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

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
    ax.scatter(p_a, p_b, p_c, s=1, color=cols)

    ax.set_tlabel('Vol. frac. A')
    ax.set_llabel('Vol. frac. B')
    ax.set_rlabel('Vol. frac. C')

    # Set axis limits
    ax.set_tlim(0, 1)
    ax.set_llim(0, 1)
    ax.set_rlim(0, 1)

    positions = ['tick1', 'tick2']
    for position in positions:
        ax.taxis.set_ticks_position(position)
        ax.laxis.set_ticks_position(position)
        ax.raxis.set_ticks_position(position)

    ax.grid()

    plt.savefig (f"signs_{chi_ab}_{chi_bc}_{chi_ac}.png", dpi=1200)
    print ("Completed heat map computation.")

