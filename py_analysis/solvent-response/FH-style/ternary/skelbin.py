import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
from scipy.optimize import root
import scipy.spatial.distance
from scipy.spatial.distance import cdist
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import mpltern
import sys
import argparse
import time
import warnings
import linecache

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}: {line}\n"

warnings.formatwarning = custom_warning_format

import argparse 
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a  memory-intensive computation.")
parser.add_argument('--chiac', metavar='chi_ac', dest='chi_ac', type=float, action='store', help='enter A-C exchange parameter.')
parser.add_argument('--chiab', metavar='chi_ab', dest='chi_ab', type=float, action='store', help='enter A-B exchange parameter.')
parser.add_argument('--chibc', metavar='chi_bc', dest='chi_bc', type=float, action='store', help='enter B-C exchange parameter.')
parser.add_argument('--mesh', metavar='mesh', dest='mesh', type=int, action='store', help='enter mesh fineness.')
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
parser.add_argument('--skelfile', dest='skelfile', type=str, action='store', help="name of file where we are going to dump the skeleton of the binodal.")
args = parser.parse_args() 



if __name__=="__main__":

    start = time.time()
    
    N = args.N
    chi_ac = args.chi_ac
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc

    va = 1
    vb = N
    vc = 1
    
    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )
    
    mesh  = args.mesh
    # phi_b_list = [np.logspace (-20, np.log10(0.9), mesh), np.linspace(0.1, 0.6, mesh), np.logspace(-15, -1, mesh)]
    phi_b_list = [np.logspace (-20, np.log10(0.9), mesh), np.hstack((np.linspace(0.6, 0.999, mesh//2), np.logspace(-20,-1,mesh//2))), np.logspace(-15, -1, mesh)]


    for idx, phi_b in enumerate(phi_b_list):

        phi_b = np.repeat (phi_b, mesh)
        phi_a = np.zeros  (phi_b.shape)

        for i in range (mesh):
            upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
            phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

        # only keep stuff which is outside the spinodal
        to_keep = stab_crit (phi_a, phi_b, chi_ab, chi_bc, chi_ac) > 0

        if np.sum (to_keep) == len(phi_a):
            print ("There is no unstable spinodal region. There will be no binodal. Moving to the next meshing of space...")
            continue

        phi_b   = phi_b [to_keep]
        phi_a   = phi_a [to_keep]

        print ("Calculating potentials...", flush=True)
        chem_pot_a = mu_a (phi_a, phi_b)
        chem_pot_b = mu_b (phi_a, phi_b)
        chem_pot_c = mu_c (phi_a, phi_b)
        print ("Calculated potentials!", flush=True)

        phis       = np.array([phi_a, phi_b, 1-phi_a-phi_b]).T
        mu         = np.array([chem_pot_a, chem_pot_b, chem_pot_c]).T 

        print ("Calculate intra-array distances...", flush=True)
        npoints = np.sum(to_keep)
        print (f"number of available points = {npoints}", flush=True)

        block_size = 500
        block_num  = np.sum (to_keep) // block_size + 1

        phi_big_block = phis[np.newaxis, :, :]
        phi_dists     = np.zeros ((npoints, npoints))

        mu_big_block  = mu[np.newaxis, :, :]
        mu_dists      = np.zeros ((npoints, npoints))
        print (f"Total number of blocks = {block_num}.")
        start_blocking = time.time()
        for i in range (block_num):
            print (f"On block number = {i}...", flush=True)
            if i < block_num - 1:
                phi_dists[i*block_size:(i+1)*block_size,:] = np.linalg.norm (phi_big_block - phis[i*block_size:(i+1)*block_size, np.newaxis,:], axis=-1)
                mu_dists [i*block_size:(i+1)*block_size,:] = np.linalg.norm (mu_big_block  - mu  [i*block_size:(i+1)*block_size, np.newaxis,:], axis=-1)
            else:
                phi_dists[i*block_size:,:] = np.linalg.norm (phi_big_block - phis[i*block_size:, np.newaxis, :], axis=-1)
                mu_dists [i*block_size:,:] = np.linalg.norm (mu_big_block  - phis[i*block_size:, np.newaxis, :], axis=-1)

        print ("Done!", flush=True)
        stop_blocking = time.time()

        print (f"Time for blocking is {stop_blocking-start_blocking} seconds.", flush=True)

        mask = phi_dists > 0.1
        print ("Calculated intra-array distances!", flush=True)

        print ("Process distances...", flush=True)
        phi_dists = np.where (mask, phi_dists, np.inf)
        mu_dists  = np.where (mask, mu_dists , np.inf)

        row_indices = np.argmin (mu_dists, axis=0)
        col_indices = np.arange (len(phi_b))
        print ("Processed!", flush=True)

        if idx == 0:
            f = open (args.skelfile, 'w')
            f.write ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
                .format ("i_1", "i_2", "dmu", "mu_a1", "mu_b1", "mu_c1", "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2")) 
        else:
            f = open (args.skelfile, 'a')


        for i in range(len(phi_b)):
            f.write  ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
                .format(col_indices[i], row_indices[i], mu_dists[row_indices[i], col_indices[i]], chem_pot_a[col_indices[i]], chem_pot_b[col_indices[i]], chem_pot_c[col_indices[i]], phis[col_indices[i],0], phis[col_indices[i],1], phis[col_indices[i],2], chem_pot_a[row_indices[i]], chem_pot_b[row_indices[i]], chem_pot_c[row_indices[i]], phis[row_indices[i],0], phis[row_indices[i],1], phis[row_indices[i],2] ) )

        f.close ()
        del phi_big_block
        del mu_big_block
        del phi_dists
        del mu_dists


    stop = time.time()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)

