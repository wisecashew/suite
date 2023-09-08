#!/home/satyend/.conda/envs/phase/bin/python

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
import multiprocessing as mp
import itertools

np.set_printoptions (threshold=sys.maxsize)

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush()


import argparse
parser = argparse.ArgumentParser(description="Create diagrams where chi_ab, chi_bc are held fixed for a plot while chi_ac is varied, serially.")
parser.add_argument('-N', metavar='N', dest='N', type=int, action='store', help='degree of polymerization of B.')
parser.add_argument('--llim', metavar='llim', dest='llim', type=float, action='store', help='lower limit of chi.')
parser.add_argument('--ulim', metavar='ulim', dest='ulim', type=float, action='store', help='upper limit of chi.')
parser.add_argument('--chunk', action='store_true', default=False, help='perform computation with chunking.')
parser.add_argument('--chunk-size', metavar='chunksize', dest='chunksize', action='store', type=int, help='Determine size of the square chunk you want from the chi mesh.', default=100)
parser.add_argument('--phimeshsize', metavar='phimeshsize', dest='phimeshsize', action='store', type=int, help='fineness of phimesh. default=100.', default=100)
parser.add_argument('--chimeshsize', metavar='chimeshsize', dest='chimeshsize', action='store', type=int, help='fineness of chimesh. default=100.', default=100)
parser.add_argument('--nproc', metavar='nproc', dest='nproc', action='store', type=int, help='enter number of processes you want to spawn.')
parser.add_argument('--chiac', metavar='chiac', dest='chiac', action='store', type=float, help='enter solvent-cosolvent mixing capacity.')
parser.add_argument('--img',   metavar='img',   dest='img',   action='store', type=str, help='enter name of image to be created.')
args = parser.parse_args() 

def stab_crit (p_a, p_b, N, c_ab, c_bc, c_ac):
    return (1/(N*p_b) + 1/(1-p_a - p_b) - 2 * c_bc) * (1/p_a + 1/(1-p_a - p_b) - 2 * c_ac) - (1/(1-p_a-p_b) + c_ab - c_bc - c_ac) ** 2


def mesh_maker (ij_tups, chunksize, chimeshsize, chi_ab, chi_bc, chi_ac, p_a, p_b):

    to_plot = []
    for elem in ij_tups:
        L1 = elem[0]*chunksize
        L2 = (elem[0]+1)*chunksize+1
        R1 = elem[1]*chunksize
        R2 = (elem[1]+1)*chunksize+1
        # print (f"L1 = {L1}, L2 = {L2}, R1 = {R1}, R2 = {R2}")
        chi_ab_s = chi_ab[L1:L2, R1:R2]
        chi_bc_s = chi_bc[L1:L2, R1:R2]
        results = stab_crit (p_a, p_b, N, chi_ab_s[:, :, np.newaxis], chi_bc_s[:, :, np.newaxis], chi_ac)
        Z = np.min (results, axis=2)
        Z = (Z<0) * 1
        to_plot.append ( (chi_ab_s, chi_bc_s, Z) )
        # if (elem[0]==0) and (elem[1]==0):
        # print (chi_ab_s)
        # print (chi_bc_s)
        # print (f"i={elem[0]}, j={elem[1]}, \nZ= {Z}")

    return to_plot

if __name__=="__main__":

    start = time.time()

    N     = args.N
    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    chi_ac = [args.chiac]# np.linspace (-10, 0, 11)

    print (f"chi_ac [0] = {chi_ac[0]}")

    chimeshsize = args.chimeshsize
    llim = args.llim
    ulim = args.ulim
    chi_ab, chi_bc = np.meshgrid ( np.linspace(llim, ulim, chimeshsize), np.linspace (llim, ulim, chimeshsize) )


    print (f"chi_ab = \n{chi_ab.shape}")
    print (f"chi_bc = \n{chi_bc.shape}")

    fig = plt.figure ()
    ax  = plt.axes   ()

    ax.axvline (x=0, c='k', ls='--')
    ax.axhline (y=0, c='k', ls='--')

    phimeshsize = args.phimeshsize
    p_a_space = np.logspace (-3, np.log10(1-0.001), phimeshsize)
    p_a = np.repeat (p_a_space, len(p_a_space))
    p_b = np.empty (p_a.shape)
    for i in range (len(p_a_space)):
        lowlim = -3 if 1-p_a_space[i] > 0.001 else -4
        uplim  = np.log10(1-p_a_space[i]-0.001)
        p_b [i*len(p_a_space):(i+1)*len(p_a_space)] = np.logspace (lowlim, uplim, len(p_a_space))

    Z = np.zeros (chi_ab.shape)

    # this is the square chunk size 
    chunksize = args.chunksize

    # max possible number of chunks, given chunk size is...
    max_chunks = (chimeshsize ** 2) // (chunksize ** 2)

    if max_chunks < args.nproc:
        nproc       = max_chunks
        cperprocess = 1
    else:
        nproc = args.nproc
        cperprocess = 1

    print (f"num process = {nproc}.")
    print (f"Chi mesh length = {chimeshsize}.")

    # check if total number of chunks can be equally distributed among all processes

    print (f"chunks assigned per process = {cperprocess}")

    ij_tups = []
    val = int(np.sqrt (nproc))
    for i in range (val):
       for j in range (val):
           ij_tups.append ([(i,j)])

    pool = mp.Pool (processes=nproc)

    results = pool.starmap (mesh_maker, zip (ij_tups, itertools.repeat(chunksize), itertools.repeat(chimeshsize), itertools.repeat(chi_ab), itertools.repeat(chi_bc), itertools.repeat(chi_ac[0]), itertools.repeat (p_a), itertools.repeat(p_b) ) )

    pool.close ()
    pool.join  ()

    for ind, indiv_proc in enumerate(results):
        for plots in indiv_proc:
            if len (plots) == 0:
                pass
            else:
                # print (f"Ind = {ind}")
                # print (f"plots[0].shape = {plots[0].shape}, plots[1].shape = {plots[1].shape}, plots[2].shape = {plots[2].shape}")
                ax.pcolormesh (plots[0], plots[1], plots[2], cmap=cm.ocean, norm=colors.Normalize(vmin=0, vmax=1), shading="auto")


    plt.savefig (args.img, dpi=1200)
    print ("Completed cononsolvency computation.")
    stop = time.time()
    print (f"Time for computation is {stop-start} seconds.")
