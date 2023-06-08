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
import sys
import argparse
import time
import os

np.set_printoptions (threshold=sys.maxsize)

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush()


import argparse
parser = argparse.ArgumentParser(description="Create diagrams where chi_ab, chi_bc are held fixed for a plot while chi_ac is varied, serially.")
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
args = parser.parse_args() 

def stab_crit (p_a, p_b, N, c_ab, c_bc, c_ac):
    return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

def get_colormesh (p_a, p_b, N, chi_ab, chi_bc, chi_ac):

    vals = stab_crit (p_a[:, np.newaxis], p_b[:, np.newaxis], N, chi_ab.ravel(), chi_bc.ravel(), chi_ac)
    Z = np.reshape ((vals<0).any(axis=0), chi_ab.shape).astype(int)

    return Z


if __name__=="__main__":

    start = time.time()

    N     = args.N
    lsize = 3
    plt.rcParams ['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    chi_ac = np.linspace (-10, 0, 11)

    chimeshsize = 25
    plot_lim = 5
    chi_ab, chi_bc = np.meshgrid ( np.linspace(-plot_lim, 0, chimeshsize), np.linspace (-plot_lim, 0, chimeshsize) )

    fig, ax = plt.subplots (figsize=(8,6), constrained_layout=True)

    phimeshsize = 25
    p_a_space = np.logspace (-3, np.log10(1-0.001), phimeshsize)
    p_a = np.repeat (p_a_space, len(p_a_space))
    p_b = np.empty (p_a.shape)
    for i in range (len(p_a_space)):
        lowlim = -3 if 1-p_a_space[i] > 0.001 else -4
        uplim  = np.log10(1-p_a_space[i]-0.001)
        p_b [i*len(p_a_space):(i+1)*len(p_a_space)] = np.logspace (lowlim, uplim, len(p_a_space))

    Z = np.zeros (chi_ab.shape)

    # broke loop method
    for i in range (chimeshsize):
        for j in range (chimeshsize):
            vals = stab_crit (p_a, p_b, N, chi_ab[i,j], chi_bc[i,j], chi_ac[0])
            print (f"(i,j) = {i,j}")
            print ((vals<0).any() * 1)
            Z[i,j] = (vals<0).any() * 1


    ax.pcolormesh     (chi_ab, chi_bc, Z, cmap=cm.ocean, norm=colors.Normalize(vmin=0, vmax=1), shading="auto")


    plt.savefig (f"cnscheck_srl", dpi=1200)
    print ("Completed cononsolvency computation.")
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.")
