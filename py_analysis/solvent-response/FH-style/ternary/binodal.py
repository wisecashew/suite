#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
import matplotlib.patheffects as pe 
from scipy.optimize import fsolve
from scipy.optimize import root
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import mpltern
import sys
import argparse
np.set_printoptions(threshold=sys.maxsize)

import argparse 
parser = argparse.ArgumentParser(description="Create a ternary spinodal diagram.")
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
    chi_ac = args.chi_ac
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc

    va = 1
    vb = 1
    vc = 1

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b, phi_c: np.log(phi_a) + 1 - phi_a - va/vb * phi_b - va/vc * phi_c + va * (phi_b**2 * chi_ab + phi_c**2 * chi_ac + phi_b * phi_c * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b, phi_c: np.log(phi_b) + 1 - phi_b - vb/va * phi_a - vb/vc * phi_c + vb * (phi_a**2 * chi_ab + phi_c**2 * chi_bc + phi_a * phi_c * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b, phi_c: np.log(phi_c) + 1 - phi_c - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    # generate a fine mesh of phi_a
    phi_c_list = np.arange (0.001, 0.999, 0.00001)
    phi_b_1 = np.zeros (phi_c_list.shape)
    phi_b_2 = np.zeros (phi_c_list.shape)

    seed_scalars = [0.01, 0.1, 0.15, 0.2]
    for idx, phi_c in enumerate(phi_c_list):
        cphi = phi_c
        def equations (phi):

            eq1 = mu_a (1-phi[0]-cphi, phi[0], cphi) - mu_a (1-phi[1]-cphi, phi[1], cphi)
            eq2 = mu_b (1-phi[0]-cphi, phi[0], cphi) - mu_b (1-phi[1]-cphi, phi[1], cphi)
            
            return [eq1, eq2]

        for seed in seed_scalars:
            sol = fsolve (equations, [(1-phi_c)*seed, (1-phi_c)*(1-seed)])
            if np.isinf(equations(sol)).any() or np.isnan(equations(sol)).any():
                pass
            elif abs (sol[0]-sol[1]) < 1e-6 or (np.abs (equations(sol)) > 1e-6).any():
                pass
            else:
                phi_b_1[idx] = sol[0]
                phi_b_2[idx] = sol[1]
                break
        

    to_keep = phi_b_1 > 0
    # Plot the points
    ax.scatter((1-phi_c_list-phi_b_1)[to_keep], phi_b_1[to_keep], phi_c_list[to_keep], s=1, c='steelblue')
    ax.scatter((1-phi_c_list-phi_b_2)[to_keep], phi_b_2[to_keep], phi_c_list[to_keep], s=1, c='steelblue')

    
    # generate a fine mesh of phi_a
    phi_a_list = np.arange (0.001, 0.999, 0.00001)
    phi_b_1 = np.zeros (phi_a_list.shape)
    phi_b_2 = np.zeros (phi_a_list.shape)

    for idx, phi_a in enumerate(phi_a_list):
        aphi = phi_a
        def equations (phi):
            eq1 = mu_b (aphi, phi[0], 1-phi[0]-aphi) - mu_b (aphi, phi[1], 1-phi[1]-aphi)
            eq2 = mu_c (aphi, phi[0], 1-phi[0]-aphi) - mu_c (aphi, phi[1], 1-phi[1]-aphi)
            
            return [eq1, eq2]

        for seed in seed_scalars:
            sol = fsolve (equations, [(1-phi_a)*seed, (1-phi_a)*(1-seed)])
            if np.isinf(equations(sol)).any() or np.isnan(equations(sol)).any():
                pass
            elif abs (sol[0]-sol[1]) < 1e-6 or (np.abs (equations(sol)) > 1e-6).any():
                pass
            else:    
                phi_b_1[idx] = sol[0]
                phi_b_2[idx] = sol[1]
                break

    to_keep = phi_b_1 > 0

    # Plot the points
    ax.scatter(phi_a_list[to_keep], phi_b_1[to_keep], (1-phi_a_list-phi_b_1)[to_keep], s=1, c='steelblue')
    ax.scatter(phi_a_list[to_keep], phi_b_2[to_keep], (1-phi_a_list-phi_b_2)[to_keep], s=1, c='steelblue')


    # generate a fine mesh of phi_a
    phi_b_list = np.arange (0.001, 0.999, 0.00001)
    phi_a_1 = np.zeros (phi_b_list.shape)
    phi_a_2 = np.zeros (phi_b_list.shape)

    for idx, phi_b in enumerate(phi_b_list):
        bphi = phi_b
        def equations (phi):
            eq1 = mu_a (1-bphi-phi[0], bphi, phi[0]) - mu_a (1-bphi-phi[1], bphi, phi[1])
            eq2 = mu_c (1-bphi-phi[0], bphi, phi[0]) - mu_c (1-bphi-phi[1], bphi, phi[1])
            
            return [eq1, eq2]

        for seed in seed_scalars:
            sol = fsolve (equations, [(1-phi_b)*seed, (1-phi_b)*(1-seed)])
            if np.isinf(equations(sol)).any() or np.isnan(equations(sol)).any():
                pass
            elif abs (sol[0]-sol[1]) < 1e-6 or (np.abs (equations(sol)) > 1e-6).any():
                pass
            else:
                # print (f"sols = {equations(sol)} , sol = {sol}, phi = {phi_a}")
                phi_a_1[idx] = sol[0]
                phi_a_2[idx] = sol[1]
                break

    to_keep = phi_a_1 > 0

    # Plot the points
    ax.scatter(phi_a_1[to_keep], phi_b_list[to_keep], (1-phi_a_1-phi_b_list)[to_keep], s=1, c='steelblue')
    ax.scatter(phi_a_2[to_keep], phi_b_list[to_keep], (1-phi_a_2-phi_b_list)[to_keep], s=1, c='steelblue')
    
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

    ax.grid()

    plt.savefig ("binodal", dpi=1200)