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

#####################

def crit_condition (N, phi_p, phi_s, chi_sc, chi_ps, chi_pc):

    phi_c = 1-phi_p-phi_s
    t1    = 1/phi_c + 1/(phi_p*N) - 2*chi_pc
    t2    = (1/(phi_c)**2 - 1/phi_s**2)*(1/(phi_c) + 1/(phi_p*N) - 2*chi_pc) + (1/(phi_c) + 1/phi_s - 2*chi_sc)/(phi_c)**2 - 2*(1/phi_c - chi_pc - chi_sc + chi_ps)/(phi_c)**2

    u1    = (1/phi_c + 1/(phi_p*N) - 2*chi_pc)/(phi_c)**2 + (1/(phi_c)**2 - 1/(phi_p**2 * N))*(1/phi_c + 1/phi_s - 2*chi_sc) - 2*(1/phi_c + chi_ps - chi_sc - chi_pc)/phi_c**2
    u2    = 1/phi_c - chi_pc - chi_sc + chi_ps

    return t1*t2 - u1*u2


#####################

def find_crit_point (N, chi_sc, chi_ps, chi_pc):

    def send_to_fsolve_r1 (phi_s):
        phi_p_upper = root_up (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (N, phi_p_upper, phi_s, chi_sc, chi_ps, chi_pc)

    def send_to_fsolve_r2 (phi_s):
        phi_p_lower = root_lo (phi_s, chi_ps, chi_pc, chi_sc)
        return crit_condition (N, phi_p_lower, phi_s, chi_sc, chi_ps, chi_pc)

    guesses = np.linspace (0, 1, 10000)
    roots_up   = np.empty ((0,2))
    roots_down = np.empty ((0,2))

    for g in guesses:
        root = fsolve (send_to_fsolve_r1, g)

        if abs(send_to_fsolve_r1(root)) < 1e-6:

            if root >= 1 or root <= 0 or np.isnan(root):
                pass
            else:
                r_up  = root_up(root, chi_ps, chi_pc, chi_sc)[0]
                r_tup = np.array([root[0], r_up])
                if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_up):
                    pass

                elif r_tup in roots_up:
                    pass

                else:
                    if len(roots_up) == 0:
                        roots_up = np.vstack ((roots_up,r_tup))
                    else:
                        similarity = (np.linalg.norm(roots_up - r_tup, axis=1) < 1e-3).any ()
                        if similarity:
                            pass
                        else:
                            roots_up = np.vstack ((roots_up,r_tup))

        else:
            pass

    for g in guesses:
        root = fsolve (send_to_fsolve_r2, g)

        if abs(send_to_fsolve_r2(root)) < 1e-6:

            if root >= 1 or root <= 0 or np.isnan(root):
                pass
            else:
                r_lo = root_lo(root, chi_ps, chi_pc, chi_sc)[0]
                r_tup = np.array([root[0], r_lo])
                if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_lo):
                    pass

                elif r_tup in roots_down:
                    pass

                else:
                    if len(roots_down) == 0:
                        roots_down = np.vstack ((roots_down,r_tup))
                    else:
                        similarity = (np.linalg.norm(roots_down - r_tup, axis=1) < 1e-3).any ()
                        if similarity:
                            pass
                        else:
                            roots_down = np.vstack ((roots_down,r_tup))

        else:
            pass

    return roots_up, roots_down




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

def perform_sweep (phi_b, mesh, chi_ab, chi_bc, chi_ac):

    phi_b = np.repeat (phi_b, mesh)
    phi_a = np.zeros  (phi_b.shape)

    print ("Populating the meshes for phi_a...", flush=True)
    for i in range (mesh):
        upper_lim = 0.999 if phi_b[i*mesh] < 0.001 else 1- phi_b[i*mesh] - 0.001
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
    
    print ("Running the distances...")
    
    phi_dists = np.linalg.norm (phis[:, np.newaxis] - phis, axis=2)
    mu_dists  = np.linalg.norm (mu[:, np.newaxis] - mu,     axis=2)

    mask      = phi_dists > 0.1
    phi_dists = np.where (mask, phi_dists, np.inf)
    mu_dists  = np.where (mask, mu_dists,  np.inf)

    print ("Pruning...", flush=True)
    closest_indices = np.argmin(mu_dists, axis=1)
    mu_conj   = mu [closest_indices]
    phi_conj  = phis[closest_indices]
    mu_dists  = mu_dists[np.arange(len(mu)), closest_indices]

    return (phis, phi_conj, mu_dists, phi_a, phi_b)

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

    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4* N * (1 - 2* phi_s * chi_sc + 2 * phi_s ** 2 * chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    (-1 + 2 * phi_s * chi_sc + N * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc:  1 / (2*N * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: 1 - 2 * phi_s * chi_sc + N * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )


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

    print (f"edges_x = {edges_x}")
    print (f"edges_y = {edges_y}")


    roots_up, roots_down = find_crit_point (N, chi_ac, chi_ab, chi_bc)
    print (f"roots_up = {roots_up}")
    print (f"roots_down = {roots_down}")
    crits = np.vstack((roots_up, roots_down))
    print (f"critical points = {crits}")

    phi_b_spin_min = 0.001
    phi_b_spin_max = 0.8

    phi_b_list = [np.linspace(phi_b_spin_min*0.9, phi_b_spin_max+(1-phi_b_spin_max)*0.1, mesh), \
                  np.linspace(phi_b_spin_min*0.7, phi_b_spin_max+(1-phi_b_spin_max)*0.3, mesh), \
                  np.linspace(phi_b_spin_min*0.5, phi_b_spin_max+(1-phi_b_spin_max)*0.5, mesh), \
                  np.linspace(phi_b_spin_min*0.1, phi_b_spin_max+(1-phi_b_spin_max)*0.9, mesh), \
                  np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh)] # , \
                  # np.hstack((np.linspace(phi_b_spin_min*0.001, phi_b_spin_max+(1-phi_b_spin_max)*0.999, mesh//2),np.linspace(np.min(crits[:,1]), np.max(crits[:,1]), mesh//2)))]

    pool    = mp.Pool (processes=len(phi_b_list))

    results = pool.starmap(perform_sweep, zip(phi_b_list, itertools.repeat(mesh), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chi_ac) ) )

    print (f"min x = {np.min(crits[:,1])}, max x = {np.max(crits[:,1])}")

    for r in results:
        plt.scatter (r[3], r[4], c='steelblue', s=0.1)
    plt.xlabel ("frac A")
    plt.ylabel ("frac B")
    plt.savefig ('testing_guess.png', dpi=1200)

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
        .format(res[2][i], res[0][i,0], res[0][i,1], res[0][i,2], res[1][i,0], res[1][i,1], res[1][i,2] ) )

    f.close ()
    stop = time.time()
    print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)



