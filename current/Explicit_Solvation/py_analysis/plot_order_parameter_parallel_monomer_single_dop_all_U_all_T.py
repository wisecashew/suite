#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np 
import re 
import matplotlib
# matplotlib.use('Agg')
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
import cmath 

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-dop', metavar='DOP'   , dest='dop', type=int , action='store', help='enter a degree of polymerization.')
parser.add_argument('-s'  , metavar='S'     , type=int  , dest='s' , action='store', help='start parsing after this index.', default=100)
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('--ortn-file', dest='of', metavar='coords', action='store', type=str, help='Name of orientation dump file to parse information.', default='coords')
parser.add_argument('--sim-style', dest='ss', action='store', type=str, help='flory or entropy.')
parser.add_argument('--norm', dest='norm', metavar='z', action='store', type=str, help='Form of normalization for fluctuations.')
parser.add_argument('--color-scheme', dest='cs', action='store', type=int, help='Type of coloring.')
parser.add_argument('--png-name', dest='pn', metavar='imagename', type=str, action='store', help='Name of image file to be generated.', default="order_param")
parser.add_argument('--show-plot', dest='sp', action ='store_true' , help='Flag to include to see plot.') 
args = parser.parse_args() 

if args.cs == 3:
    divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0.0, vmax=1.0 ) 
elif args.cs == 0:
    if args.ss == 'flory-huggins':
        divnorm = matplotlib.colors.SymLogNorm ( 0.0001, vmin=-0.2, vmax=0.1 ) # this is for flory-huggins  
    elif args.ss == 'entropy':
        divnorm = matplotlib.colors.SymLogNorm ( 0.02, vmin=-0.2, vmax=0.1 ) # this is for entropy 
    else:
        print ("Bad simstyle.", flush=True)
        exit()

if __name__=="__main__":

    start = time.time() 
    #######################################
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes  () 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params (axis='x', labelsize=16)
    ax.tick_params (axis='y', labelsize=16)

    U_list    = aux.dir2U ( os.listdir (".") ) 
    # U_list    = ["U1", "U2", "U4", "U5", "U10"]
    PLOT_DICT = {} 
    
    starting_index = args.s
    ortn_file      = args.of
    dop            = args.dop
    show_plot_bool = args.sp
    i    = 0
    nproc = args.nproc
    max_op = 0 
    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 ) 
    pool_list = [ pool1 ]

    for U in U_list:
        print ("Inside U = "+str(U)+"...", flush=True)
        temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(dop) ) )
        # temperatures = [100.0]
        # get the num_list for each temperature 
        master_temp_list = [] 
        master_num_list  = [] 
        ord_par1_dict     = {} 
        ord_par2_dict     = {} 
        ntraj_dict       = {} 

        # define stores 
        ord_par1_mean     = [] 
        ord_par1_std      = [] 

        for T in temperatures:
            num_list         = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend  ( num_list )
            master_temp_list.extend ( [T]* len(num_list) )
            ntraj_dict[T]    = len (num_list)
            ord_par1_dict [T] = []

        # start multiprocessing... keeping in mind that each node only has 96 cores
        # start splitting up master_num_list and master_temp_list 

        # define a shitty dict
        # shitty_dict = {0:0, 1:0, 2:1}
        idx_range = len(master_num_list)//nproc + 1
        for u_idx in range(idx_range):
            if u_idx == idx_range-1:
                results = pool_list [ 0 ].starmap ( aux.obtain_monomer_alignment, zip( itertools.repeat(U), itertools.repeat(dop), master_temp_list[u_idx*nproc:], itertools.repeat(ortn_file), master_num_list[u_idx*nproc:], itertools.repeat(starting_index), itertools.repeat(args.norm) ) ) 
            else:
                results = pool_list [ 0 ].starmap ( aux.obtain_monomer_alignment, zip( itertools.repeat(U), itertools.repeat(dop), master_temp_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(ortn_file), master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(starting_index), itertools.repeat(args.norm) ) ) 
                
            print ("Pool has been closed. This pool has {} threads.".format ( len(results ) ), flush=True )

            for k in range ( len (master_temp_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
                ord_par1_dict[ master_temp_list[u_idx*nproc + k] ].append ( results[k] )

        for T in np.unique ( master_temp_list ):
            ord_par1_mean.append  ( np.mean  ( ord_par1_dict[T] ) )
            ord_par1_std.append   ( np.std   ( ord_par1_dict[T] )/ np.sqrt(master_temp_list.count(T)) )
        if max_op < np.max(ord_par1_mean):
            max_op = np.max(ord_par1_mean) 
        
        # print (ord_par1_mean)
        PLOT_DICT [U] = (np.asarray (ord_par1_mean), np.asarray (ord_par1_std) ) 
    pool1.close()
    pool1.join () 

    i = 0 
    print ("max = ",max_op)
    max_op = 1
    for U in U_list: 
        chi_a = aux.get_chi_cosolvent ( str(U)+"/geom_and_esurf.txt")[args.cs]
        rgba_color = cm.PiYG (divnorm(chi_a))
        ax.errorbar ( temperatures, PLOT_DICT[U][0]/max_op, yerr=PLOT_DICT[U][1], fmt='none', linestyle='-', elinewidth=1, capsize=2, linewidth=1, color=rgba_color, label='_nolegend_')
        ax.plot     ( temperatures, PLOT_DICT[U][0]/max_op, marker='o', markeredgecolor='k', linestyle='-', linewidth=2, c=rgba_color, label='_nolegend_', markersize=10)
        i += 1

    ##############################################################

    my_cmap = cm.PiYG
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=-0.2, vmax=0.1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1] )
    cbar.ax.tick_params (labelsize=14)
    cbar.set_ticklabels( ["-$10^{-1}$", "-$10^{-2}$","-$10^{-3}$",  0, "$10^{-3}$", "$10^{-2}$", "$10^{-1}$" ] )
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    ax.set_ylim((-0.02, 0.42))
    ax.set_xscale('log')
    ax.minorticks_on()
    ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
    ax.set_yticks (np.arange (0, 0.45, 0.1) )
    plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.1f}'))
    # ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    plt.savefig   ( args.pn + ".png", dpi=1000 )

    ##################################
    stop = time.time() 
    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

