#!/home/satyend/.conda/envs/phase/bin/python

import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import numpy as np 
import re 
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import os
import aux
import time
import sys
import multiprocessing
import itertools

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

'''
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

import argparse
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number (not index or line number in file).', default=2000)
parser.add_argument('--u', metavar='u', dest='u', action='store', nargs='+', type=str, help='Enter force-fields you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('--coords', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--png-name', dest='pn', metavar='imagename', action='store', type=str, help='Name of image file', default='rg_plot')
args = parser.parse_args() 

divnorm = matplotlib.colors.LogNorm ( vmin=0.01, vmax=100 ) 


def get_starting_ind ( U, frac, num, dumpfile, s):
    filename = U + "/" + frac + "/" + dumpfile + "_" + str(num) + ".mc"
    df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    L = len(df["energy"])
    return int(df["time_step"].values[L-s])


def infiltrate_coords_get_rg ( U, frac, num, coords_files, starting_index ):

    dop = 32
    filename    = U + "/" + frac + "/"+ coords_files + "_" + str(num) + ".mc" 
    edge        = aux.edge_length (dop)
    master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
    rg          = aux.get_Rg(master_dict, edge, edge, edge) 
    return rg



if __name__ == "__main__":

    start = time.time()

    ##################################

    def mysorter_f (x):
        new_str = ""
        for m in x:
            if m.isdigit():
                new_str = new_str+m
        return float(new_str)

    u_list = args.u # aux.dir2u ( os.listdir (".") )
    PLOT_DICT = {} 
    coords_files   = args.c
    starting_index = args.s

    ######

    fig = plt.figure( figsize=(5,5), constrained_layout=True )
    ax  = plt.axes  ()
    ax.tick_params  ( direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params  ( axis='x', labelsize=8)
    ax.tick_params  ( axis='y', labelsize=8)
    i = 0

    rg_max = 1

    # instantiating pool
    nproc = args.nproc
    pool1 = multiprocessing.Pool ( processes=nproc ) # len(num_list))

    pool_list = [pool1]


    for u in u_list:
        print (f"Inside u = {u}...", flush=True)
        frac_list = list (os.listdir(u))
        frac_list.sort()
        xu_list   = [u]*len(frac_list)
        rg_mean   = []
        rg_std    = []

        # get num_list for each temperature
        master_u_list     = []
        master_num_list   = []
        master_frac_list  = []
        master_index_list = []
        rg_dict    = {}
        ntraj_dict = {}

        for frac in frac_list:
            num_list = list ( np.unique ( aux.dir2nsim ( os.listdir (str(u) + "/" + frac) ) ) )
            master_num_list.extend ( num_list )
            master_frac_list.extend ([frac]*len(num_list))
            for num in num_list:
                master_index_list.append (get_starting_ind (u, frac, num, "energydump", starting_index) )
            master_u_list.extend ( [u]*len( num_list ) )
            ntraj_dict[frac] = len ( num_list )
            rg_dict[frac] = []

        # start multiprocessing... keeping in mind that each node only has 96 cores
        # start splitting up master_num_list and master_temp_list
        print (f"Preparing a pool of {nproc} threads to be launched repeatedly.", flush=True)

        idx_range = len (master_num_list)//nproc + 1
        for u_idx in range (idx_range):
            if u_idx == idx_range-1:
                results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( master_u_list[u_idx*nproc:], master_frac_list[u_idx*nproc:],  master_num_list [u_idx*nproc:], itertools.repeat (coords_files), master_index_list[u_idx*nproc:] ) )
            else:
                results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( master_u_list[(u_idx)*nproc:(u_idx+1)*nproc], master_frac_list[u_idx*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(coords_files), master_index_list[u_idx*nproc:(u_idx+1)*nproc] ) )

            print ("Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

            for k in range( len( master_frac_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
                rg_dict[master_frac_list[u_idx*nproc + k]].append( results[k] )

        # sorted_U_list = list(np.unique (master_U_list))
        # sorted_U_list.sort(key=mysorter_f)

        for f in frac_list:
            rg_mean.append( np.mean ( rg_dict[f] ) )
            rg_std.append ( np.std  ( rg_dict[f] ) / np.sqrt(master_frac_list.count(f) ) )

        PLOT_DICT [u] = (np.asarray(rg_mean), np.asarray(rg_std))


    pool_list[0].close ()
    pool_list[0].join  ()

    print ("rg_max = ", rg_max)

    toplot_frac = [float(elem[-3:]) for elem in frac_list]
    print (u_list)
    for idx,u in enumerate(u_list):
        chi = aux.get_chi_sc (u+"/"+frac_list[idx]+"/geom_and_esurf_"+frac_list[idx][-3:]+".txt")
        print (chi)
        ax.errorbar ( toplot_frac, PLOT_DICT[u][0]/rg_max, yerr= PLOT_DICT[u][1]/rg_max, linewidth=1, capsize=2, color='k', fmt='none', label='_nolegend_')
        ax.plot     ( toplot_frac, PLOT_DICT[u][0]/rg_max, marker='o', markeredgecolor='k', linestyle='-', linewidth=3/2, label=f"$\\chi _{{sc}}$ = {chi}", markersize=4 ) 
        i += 1

    print (f"frac_list = {toplot_frac}")


    #########
    # plt.rcParams['text.usetex'] = True
    ax.minorticks_on()
    ax.legend (loc="upper right")
    ax.set_xlim((-0.05, 1.05))
    # ax.set_ylim(0.4, 1.5)
    ax.set_xticks (np.linspace (0, 1, 6))
    # ax.set_yticks (np.arange(0.4, 1.7, 0.2))
    # ax.set_yticklabels (ax.get_yticks(), weight='bold')
    ax.set_xticklabels (ax.get_xticks())
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))

    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    plt.savefig   ( args.pn, bbox_inches='tight', dpi=1200 )

    ##################################

    stop = time.time()
    print ("Run time for is {:.2f} seconds.".format(stop-start), flush=True)

