import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patheffects as pe
from scipy.optimize import fsolve
from scipy.optimize import root
import matplotlib.ticker as tck
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
import sys
import argparse
import time
import multiprocessing as mp
import os
import copy
import itertools
from matplotlib.colors import LinearSegmentedColormap

np.set_printoptions (threshold=sys.maxsize)

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush()


import argparse
parser = argparse.ArgumentParser(description="Create diagrams where chi_ab, chi_bc are held fixed for a plot while chi_ac is varied, in parallel.")
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
parser.add_argument('--nproc', metavar='nproc', dest='nproc', type=int, action='store', help='number of processes to spawn.')
args = parser.parse_args() 

def stab_crit (p_a, p_b, N, c_ab, c_bc, c_ac):
    return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2

def get_colormesh (p_a, p_b, N, chi_ab, chi_bc, chi_ac):
    vals = stab_crit (p_a[:, np.newaxis], p_b[:, np.newaxis], N, chi_ab.ravel(), chi_bc.ravel(), chi_ac)
    Z = np.reshape ((vals<0).any(axis=0), chi_ab.shape).astype(int)
    return Z

def sector_processor (x, y, p_a, p_b, N, chi_ab, chi_bc, chi_ac):

    sector_chi_ab = chi_ab[x[0]:x[1], y[0]:y[1]]
    sector_chi_bc = chi_bc[x[0]:x[1], y[0]:y[1]]
    Z = get_colormesh (p_a, p_b, N, sector_chi_ab, sector_chi_bc, chi_ac)

    return [sector_chi_ab, sector_chi_bc, Z]

if __name__=="__main__":

    start = time.time()
    cols = ['white', 'coral']

    # Create a colormap with two colors
    cmap = LinearSegmentedColormap.from_list('white_to_coral', cols, N=256)

    N     = args.N
    nproc = args.nproc
    lsize = 3
    # plt.rcParams ['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    chi_ac = [-10, -7, -4, 0] # np.linspace (0, 100, 11)

    chimeshsize = 1000
    plot_lim = 15
    chi_ab, chi_bc = np.meshgrid ( np.linspace(-plot_lim, 0, chimeshsize), np.linspace (-plot_lim, 0, chimeshsize) )


    phimeshsize = 100
    p_a_space = np.logspace (-3, np.log10(1-0.001), phimeshsize)
    p_a = np.repeat (p_a_space, len(p_a_space))
    p_b = np.empty (p_a.shape)
    for i in range (len(p_a_space)):
        lowlim = -3 if 1-p_a_space[i] > 0.001 else -4
        uplim  = np.log10(1-p_a_space[i]-0.001)
        p_b [i*len(p_a_space):(i+1)*len(p_a_space)] = np.logspace (lowlim, uplim, len(p_a_space))


    # i am going to chunk up chis into 96 blocks of size nchunk
    nchunk = int(np.floor(np.sqrt (chimeshsize*chimeshsize/nproc))) + 1
    
    x_lims = list(range(0,chimeshsize,nchunk))
    x_lims.append(chimeshsize)
    y_lims = copy.copy(x_lims)
    
    x_tups = []
    y_tups = []

    for i in range(len(x_lims)-1):
        for j in range(len(y_lims)-1):
            x_tups.append ((x_lims[i],x_lims[i+1]))
            y_tups.append ((y_lims[j],y_lims[j+1]))


    # print (f"chi_ab = \n{chi_ab}")

    spawns = len(x_tups)


    pool = mp.Pool (processes=spawns)

    for idx, chiac in enumerate(chi_ac):
        fig2  = plt.figure (num=idx, figsize=(2.0,2.0))
        ax2   = plt.axes   ( )
        ax2.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')

        # ax2.axvline (x=0, c='k', linestyle='--')
        # ax2.axhline (y=0, c='k', linestyle='--')
        ax2.set_xticks([-15, -10, -5, 0])
        ax2.set_yticks([-15, -10, -5, 0])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())

        print (f"num = {idx}, chi_ac = {chiac}")
        results = pool.starmap (sector_processor, zip(x_tups, y_tups, itertools.repeat(p_a), itertools.repeat(p_b), itertools.repeat(N), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chiac) ) )
        for things in results:
            if len(things) == 0:
                pass
            else:
                ax2.pcolormesh (things[0], things[1], things[2], cmap=cmap, norm=colors.Normalize(vmin=0, vmax=1), shading="auto")
        ax2.minorticks_on()
        fig2.savefig (f"cnscheck_chiac"+str(chiac), dpi=1200, bbox_inches="tight")

    pool.close ()
    pool.join  ()

    # print (f"len(results) = {len(results)}")

    print ("Completed cononsolvency computation.")
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.")
