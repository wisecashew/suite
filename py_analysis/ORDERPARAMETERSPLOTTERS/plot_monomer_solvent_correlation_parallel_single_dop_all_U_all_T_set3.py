#!/usr/licensed/anaconda3/2020.7/bin/python

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
import time 
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux 
import multiprocessing 
import itertools
import cmath 
from pathlib import Path

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

import argparse 
parser = argparse.ArgumentParser(description="Plot monomer-solvent correlation.")
parser.add_argument('-dop', metavar='DOP'   , dest='dop', type=int , action='store', help='enter a degree of polymerization.')
parser.add_argument('-s'  , metavar='S'     , type=int  , dest='s' , action='store', help='start parsing after this index.', default=100)
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('--ortn-file', dest='of', metavar='orientation', action='store', type=str, help='Name of orientation dump file to parse information.', default='orientation')
parser.add_argument('--png-name', dest='pn', metavar='imagename', type=str, action='store', help='Name of image file to be generated.', default="order_param")

args = parser.parse_args() 

if __name__=="__main__":

    start = time.time()

    #######################################
    fig = plt.figure   ( figsize=(4/1.6,3/1.6), constrained_layout=True )
    ax  = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    norm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 ) # this is for entropy 
    font = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size':11}

    # define the colormap for the background
    # Define the gradient colors
    color1 = np.array([131, 159, 192]) / 255.0  # #839FC0 in RGB
    color2 = np.array([241, 156, 118]) / 255.0   # #ED8151 in RGB
    cmap = colors.LinearSegmentedColormap.from_list('custom', [color1, "white", color2])
    cnorm = colors.TwoSlopeNorm (vcenter=0.5, vmin=0.3, vmax=0.8)
    y    = np.array ([0,1])

    df = pd.read_csv ("INTEGRATED-FLORY-EXPONENT-TYPE2.csv", sep='|', engine='python', skiprows=1, names=["U", "T", "nu_mean", "nu_err"])
    df = df.loc[df["U"] == "U4"]
    temperatures = df["T"].values

    x_old = temperatures
    y_old = df["nu_mean"].values/2

    x_pred = np.logspace(-2, 2, 10000)
    y_pred = np.interp  (x_pred, x_old, y_old)

    col_dict = dict ()

    for i in range (len(x_pred)):
        col_dict [ x_pred [i] ] = y_pred [i]


    X, Y = np.meshgrid (x_pred, y)
    Z    = np.zeros ((2, len(x_pred)))

    for i in range (len(x_pred)):
        for j in range (2):
            Z[j, i] = col_dict [x_pred[i]]

    ax.pcolormesh (X, Y, Z, cmap=cmap, norm=cnorm, shading="auto")



    U_list    = ["U4"]







    PLOT_DICT = dict ()
    how_much_to_skip = args.s
    ortn_file      = args.of
    dop            = args.dop

    i    = 0
    nproc = args.nproc
    max_op = 0

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=args.nproc ) 
    pool_list = [ pool1 ]

    for U in U_list:
        print ("Inside U = "+str(U)+"...", flush=True)
        # temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(dop) ) )
        temperatures = [0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
        print (temperatures)

        # temperatures = [0.01]
        # get the num_list for each temperature
        master_temp_list  = []
        master_num_list   = []
        master_index_list = []
        ord_par1_dict     = {}
        ord_par2_dict     = {}
        ntraj_dict        = {}

        # define stores
        ord_par1_mean     = []
        ord_par1_std      = []

        for T in temperatures:
            num_list         = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend  ( num_list )
            master_temp_list.extend ( [T]* len(num_list) )
            ntraj_dict[T]    = len (num_list)
            # for num in num_list:
            #    master_index_list.append (get_starting_ind (U, T, num, dop, "energydump") )
            ord_par1_dict [T] = []

        # start multiprocessing... keeping in mind that each node only has 96 cores
        # start splitting up master_num_list and master_temp_list 

        idx_range = len(master_num_list)//nproc + 1
        for u_idx in range(idx_range):
            if u_idx == idx_range-1:
                results = pool_list [ 0 ].starmap ( aux.obtain_correlation, zip( itertools.repeat(U), itertools.repeat(dop), master_temp_list[u_idx*nproc:], itertools.repeat(ortn_file), master_num_list[u_idx*nproc:], itertools.repeat(how_much_to_skip) ) )
            else:
                results = pool_list [ 0 ].starmap ( aux.obtain_correlation, zip( itertools.repeat(U), itertools.repeat(dop), master_temp_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(ortn_file), master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(how_much_to_skip) ) )
                
            print ("Pool has been closed. This pool has {} threads.".format ( len(results ) ), flush=True )
            # print ("T = ", master_temp_list[u_idx*nproc])
            # print ("results = ",results)
            # print ("len(master_temp_list) = {}".format (len(master_temp_list)))
            for k in range ( len (master_temp_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
                # print ("k = {}".format(k), end=", ")
                ord_par1_dict[ master_temp_list[u_idx*nproc + k] ].append ( results[k] )

        for T in np.unique ( master_temp_list ):
            ord_par1_mean.append  ( np.mean  ( ord_par1_dict[T] ) )
            ord_par1_std.append   ( np.std   ( ord_par1_dict[T] ) / np.sqrt(master_temp_list.count(T) ) )
        if max_op < np.max(ord_par1_mean):
            max_op = np.max(ord_par1_mean)

        # print (ord_par1_mean)
        PLOT_DICT [U] = (np.asarray (ord_par1_mean), np.asarray (ord_par1_std) ) 
    pool1.close()
    pool1.join () 

    i = 0
    max_op = 25*2+(dop-2)*24
    for U in U_list: 
        chi_a = 0.1
        rgba_color = '#B9B41F'  # cm.PiYG (norm(chi_a))
        # df = pd.read_csv ("INTEGRATED-FLORY-EXPONENT-TYPE2.csv", sep='|')
        # nu = df[df["U"]==U]
        # nu = df.loc[df["T"].isin(temperatures)]
        ax.errorbar ( temperatures, PLOT_DICT[U][0]/max_op, yerr=PLOT_DICT[U][1]/max_op, fmt='none', linestyle='-', elinewidth=1, capsize=2, linewidth=1, color=rgba_color, label='_nolegend_', zorder=10, clip_on=False)
        ax.plot     ( temperatures, PLOT_DICT[U][0]/max_op, marker='o', markeredgecolor='k', linestyle='-', linewidth=2, c=rgba_color, label='_nolegend_', markersize=10, zorder=10, clip_on=False)
        i += 1

    ##############################################################

    # ax.minorticks_on() 
    ax.set_xscale('log')
    yticks = np.arange(0.0, 1.2, 0.2)
    ax.set_yticks ( yticks )
    ax.set_yticklabels (ax.get_yticks(), fontdict=font)
    ax.set_ylim   ( 0.0, 1.0 )
    ax.set_xlim   ( 0.01, 100 )
    ax.set_xticks (np.logspace(-2, 2, 5))
    ax.set_xticklabels (["0.01", "0.1", "1.0", "10.0", "100.0"], fontdict=font)
    ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}') )
    ax.set_aspect ('auto')
    plt.savefig   ( args.pn + ".png", dpi=1200, bbox_inches='tight')

    ##################################
    stop = time.time() 
    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

