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
import sys
import argparse
import time
import warnings
import linecache
import multiprocessing as mp
import itertools

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}: {line}"

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


def process_sweep(args):
    try:
        return perform_sweep(*args)
    except Exception as e:
        traceback.print_exc()
        return []



def perform_sweep (phi_b, mesh, chi_ab, chi_bc, chi_ac):

    phi_b = np.repeat (phi_b, mesh)
    phi_a = np.zeros  (phi_b.shape)

    for i in range (mesh):
        upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1-phi_b[i*mesh] - 0.001
        phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, upper_lim, mesh)

    # only keep stuff which is outside the spinodal
    to_keep = stab_crit (phi_a, phi_b, chi_ab, chi_bc, chi_ac) > 0

    if np.sum (to_keep) == len(phi_a):
        return []

    phi_b   = phi_b [to_keep]
    phi_a   = phi_a [to_keep]

    chem_pot_a = mu_a (phi_a, phi_b)
    chem_pot_b = mu_b (phi_a, phi_b)
    chem_pot_c = mu_c (phi_a, phi_b)

    phis       = np.array([phi_a, phi_b, 1-phi_a-phi_b]).T
    mu         = np.array([chem_pot_a, chem_pot_b, chem_pot_c]).T 

    npoints = np.sum(to_keep)

    block_size = 500
    block_num  = np.sum (to_keep) // block_size + 1

    phi_big_block = phis[np.newaxis, :, :]
    phi_dists     = np.zeros ((npoints, npoints))

    mu_big_block  = mu[np.newaxis, :, :]
    mu_dists      = np.zeros ((npoints, npoints))

    for i in range (block_num):
    
        if i < block_num - 1:
            phi_dists[i*block_size:(i+1)*block_size,:] = np.linalg.norm (phi_big_block - phis[i*block_size:(i+1)*block_size, np.newaxis,:], axis=-1)
            mu_dists [i*block_size:(i+1)*block_size,:] = np.linalg.norm (mu_big_block  - mu  [i*block_size:(i+1)*block_size, np.newaxis,:], axis=-1)
        else:
            phi_dists[i*block_size:,:] = np.linalg.norm (phi_big_block - phis[i*block_size:, np.newaxis, :], axis=-1)
            mu_dists [i*block_size:,:] = np.linalg.norm (mu_big_block  - phis[i*block_size:, np.newaxis, :], axis=-1)

    mask = phi_dists > 0.1
    phi_dists = np.where (mask, phi_dists, np.inf)
    mu_dists  = np.where (mask, mu_dists , np.inf)

    row_indices = np.argmin (mu_dists, axis=0)
    col_indices = np.arange (len(phi_b))

    del phi_big_block
    del mu_big_block
    del phi_dists

    return (col_indices, row_indices, mu_dists, chem_pot_a, chem_pot_b, chem_pot_c, phis)


if __name__=="__main__":

    start = time.time()

    N      = args.N
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
    phi_b_list = [np.logspace (-20, np.log10(0.999), mesh), np.linspace(0.1, 0.6, mesh), np.linspace(0.6, 0.999, mesh), np.logspace(-15, -1, mesh)]

    pool    = mp.Pool (processes=len(phi_b_list) )

    results = pool.starmap(perform_sweep, zip(phi_b_list, itertools.repeat(mesh), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chi_ac) ) )

    pool.close ()
    pool.join  ()

    f = open (args.skelfile, 'w')
    f.write ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
    .format ("i_1", "i_2", "dmu", "mu_a1", "mu_b1", "mu_c1", "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2")) 

    for idx,res in enumerate(results):
        if len(res) == 0:
            print (f"No critical condition in process {idx}.")
            continue
        for i in range(len (phi_b_list[idx]) ):
            f.write  ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
        .format(res[0][i], res[1][i], res[2][res[1][i], res[0][i]], res[3][res[0][i]], res[4][res[0][i]], res[5][res[0][i]], res[6][res[0][i],0], res[6][res[0][i],1], res[6][res[0][i],2], res[3][res[1][i]], res[4][res[1][i]], res[5][res[1][i]], res[6][res[1][i],0], res[6][res[1][i],1], res[6][res[1][i],2] ) )

    f.close ()
    stop = time.time()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)

