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
import numba as nb
from numba import jit

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

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

######
def process_sweep(args):
    try:
        return perform_sweep(*args)
    except Exception as e:
        traceback.print_exc()
        return []
######

def get_limits (N):

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

    to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0) * (r1+phi_a<=1)
    r1 = r1[to_keep_1]

    to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0) * (r2+phi_a<=1)
    r2 = r2[to_keep_2]

    if len(r1) == 0 or len(r2) == 0:
        return None, None

    else:
        curve_x = np.hstack((phi_a[to_keep_1], phi_a[to_keep_2]))
        curve_y = np.hstack((r1, r2))
        edges_x = [np.min(curve_x), np.max (curve_x)]
        edges_y = [np.min(curve_y), np.max (curve_y)]
        return edges_x, edges_y

#########################################

from numba import njit, prange

@njit(parallel=True)
def perform_sweep(phi_b, N, mesh, chi_ab, chi_bc, chi_ac):
    
    va = 1; vb = N; vc = 1;

    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )
    phi_b = np.repeat(phi_b, mesh)
    phi_a = np.zeros_like(phi_b)

    phi_b_mask = phi_b < 0.001
    upper_lim = np.where(phi_b_mask, 0.999, 1 - phi_b - 0.001)
    
    for i in prange(len(phi_b) // mesh):
        start_idx = i * mesh
        end_idx = (i + 1) * mesh
        phi_a[start_idx:end_idx] = np.linspace(0.001, upper_lim[start_idx], mesh)

    to_keep = stab_crit(phi_a, phi_b, chi_ab, chi_bc, chi_ac) > 0

    if np.all(to_keep):
        return []

    phi_a = phi_a[to_keep]
    phi_b = phi_b[to_keep]

    chem_pot_a = mu_a(phi_a, phi_b)
    chem_pot_b = mu_b(phi_a, phi_b)
    chem_pot_c = mu_c(phi_a, phi_b)

    phis = np.column_stack((phi_a, phi_b, 1 - phi_a - phi_b))
    mu = np.column_stack((chem_pot_a, chem_pot_b, chem_pot_c))

    npoints = len(phi_b)

    phi_dists = np.linalg.norm(phis[:, np.newaxis, :] - phis[:, np.newaxis, :], axis=-1)
    mu_dists = np.linalg.norm(mu[:, np.newaxis, :] - mu[:, np.newaxis, :], axis=-1)

    mask = phi_dists > 0.1
    phi_dists = np.where(mask, phi_dists, np.inf)
    mu_dists = np.where(mask, mu_dists, np.inf)

    row_indices = np.argmin(mu_dists, axis=0)
    col_indices = np.arange(len(phi_b))

    return col_indices, row_indices, mu_dists, chem_pot_a, chem_pot_b, chem_pot_c, phis

#########################################

if __name__=="__main__":

    start = time.time()

    N      = args.N
    chi_ac = args.chi_ac
    chi_ab = args.chi_ab
    chi_bc = args.chi_bc

    va = 1
    vb = N
    vc = 1

    @njit(parallel=True)
    def stab_crit (p_a, p_b, c_ab, c_bc, c_ac):
        return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )

    mesh  = args.mesh

    edges_x, edges_y = get_limits (N)

    if edges_x == None or edges_y == None:
        print ("There is no critical region in this regime. Terminating calculation.")
        exit  ()

    phi_b_spin_min = edges_y[0]
    phi_b_spin_max = edges_y[1]

    print (phi_b_spin_min)
    print (phi_b_spin_max)
    print (phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1)
    print (phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3)
    print (phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5)
    print (phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9)
    print (phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999)

    phi_b_list = [np.linspace(phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1, mesh), \
                  np.linspace(phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3, mesh), \
                  np.linspace(phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5, mesh), \
                  np.linspace(phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9, mesh), \
                  np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh), \
                  np.hstack((np.linspace(0.001, 0.5, mesh), np.linspace(0.7, 0.9, mesh)))]

    pool    = mp.Pool (processes=len(phi_b_list))

    results = pool.starmap(perform_sweep, zip(phi_b_list, itertools.repeat(N), itertools.repeat(mesh), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chi_ac) ) )

    pool.close ()
    pool.join  ()

    f = open (args.skelfile, 'w')
    f.write ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} \n"\
    .format ("dmu", "phi_a1", "phi_b1", "phi_c1", "phi_a2", "phi_b2", "phi_c2")) 

    for idx,res in enumerate(results):
        if len(res) == 0:
            print (f"No critical condition in process {idx}.")
            continue
        for i in range(len (phi_b_list[idx]) ):
            f.write  ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} \n"\
        .format(res[2][res[1][i], res[0][i]], res[6][res[0][i],0], res[6][res[0][i],1], res[6][res[0][i],2], res[6][res[1][i],0], res[6][res[1][i],1], res[6][res[1][i],2] ) )

    f.close ()
    stop = time.time()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)

