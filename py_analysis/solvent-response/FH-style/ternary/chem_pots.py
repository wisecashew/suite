#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
import matplotlib.patheffects as pe 
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.spatial.distance import cdist 
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import mpltern
import sys
import argparse
import time
# np.set_printoptions(threshold=sys.maxsize)

import argparse 
parser = argparse.ArgumentParser(description="Create a ternary spinodal diagram.")
parser.add_argument('--chiac', metavar='chi_ac', dest='chi_ac', type=float, action='store', help='enter A-C exchange parameter.')
parser.add_argument('--chiab', metavar='chi_ab', dest='chi_ab', type=float, action='store', help='enter A-B exchange parameter.')
parser.add_argument('--chibc', metavar='chi_bc', dest='chi_bc', type=float, action='store', help='enter B-C exchange parameter.')
parser.add_argument('--mesh', metavar='mesh', dest='mesh', type=int, action='store', help='enter mesh fineness.')
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
parser.add_argument('--dumpfile', dest='dumpfile', type=str, action='store', help="name of dumpfile.")
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

    # FIND PHI_B GIVEN PHI_A
    mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - va/vb * phi_b - va/vc * (1-phi_a-phi_b) + va * (phi_b**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_ac + phi_b * (1-phi_a-phi_b) * (chi_ab + chi_ac - chi_bc) ) 
    mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vb/va * phi_a - vb/vc * (1-phi_a-phi_b) + vb * (phi_a**2 * chi_ab + (1-phi_a-phi_b)**2 * chi_bc + phi_a * (1-phi_a-phi_b) * (chi_ab + chi_bc - chi_ac) )
    mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/va * phi_a - vc/vb * phi_b + vc * (phi_a**2 * chi_ac + phi_b**2 * chi_bc + phi_a * phi_b * (chi_ac + chi_bc - chi_ab) )
    

    # generate a fine mesh 
    mesh  = args.mesh
    phi_b = np.linspace (0.001, 0.999, mesh)
    phi_b = np.repeat (phi_b, mesh)
    print (phi_b.shape)
    phi_a = np.zeros  (phi_b.shape)
    for i in range (mesh):
        phi_a[i*mesh:(i+1)*mesh] = np.linspace (0.001, 1-phi_b[i*mesh]-0.001, mesh)

    print ("Calculating potentials...")
    chem_pot_a = mu_a (phi_a, phi_b)
    chem_pot_b = mu_b (phi_a, phi_b)
    chem_pot_c = mu_c (phi_a, phi_b)
    print ("Calculated potentials!")

    phis       = np.array([phi_a, phi_b, 1-phi_a-phi_b]).T
    chem_pot   = np.array([chem_pot_a, chem_pot_b, chem_pot_c]).T 


    phi_distance      = cdist (phis, phis)
    chem_pot_distance = cdist (chem_pot, chem_pot)
    mask = phi_distance > 0.1

    phi_distance = np.where (mask, phi_distance, np.inf)
    chem_pot_distance = np.where (mask, chem_pot_distance, np.inf)

    # print (f"phi_distance.shape = {phi_distance.shape}")
    # print (f"chem_pot_distance.shape = {chem_pot_distance.shape}")


    # now, only keep the minimum distances out of all the ones that passed that filter
    # min_chem_pot_dist     = np.amin   (chem_pot_distance, axis=0)
    row_indices = np.argmin (chem_pot_distance, axis=0)
    # print (f"row_indices = {row_indices}")
    col_indices = np.arange (mesh*mesh)


    print ("Create comparing_mu.txt...")

    f = open ("comparing_mu.txt", 'w')
    f.write ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
        .format ("i_1", "i_2", "dmu", "mu_a1", "mu_b1", "mu_c1", "phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2")) # i_1 i_2 dmu phi_a1 phi_a2 m_a1 m_a2 phi_b1 phi_b2 m_b1 m_b2 phi_c1 phi_c2 m_c1 m_c2\n")
    for i in range(mesh*mesh):
        # print (col_indices[i], row_indices[i]) #chem_pot_distance[row_indices[i], col_indices[i]], phis[col_indices[i],0], phis[row_indices[i],0], phis[col_indices[i],1], phis[row_indices[i],1], phis[col_indices[i],2], phis[row_indices[i],2])
        f.write  ("{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n"\
            .format(col_indices[i], row_indices[i], chem_pot_distance[row_indices[i], col_indices[i]], chem_pot_a[col_indices[i]], chem_pot_b[col_indices[i]], chem_pot_c[col_indices[i]], phis[col_indices[i],0], phis[col_indices[i],1], phis[col_indices[i],2], chem_pot_a[row_indices[i]], chem_pot_b[row_indices[i]], chem_pot_c[row_indices[i]], phis[row_indices[i],0], phis[row_indices[i],1], phis[row_indices[i],2] ) )

    f.close ()
    print ("Created!")
    stop = time.time() 
    print (f"Time elapsed is {stop-start} seconds.")
    
    print ("Completed standard ternary plots.")