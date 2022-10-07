#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np 
import re 
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
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
parser.add_argument('--ortn-file', dest='of', metavar='orientation', action='store', type=str, help='Name of orientation dump file to parse information.', default='orientation')
parser.add_argument('--ord-param', dest='op', metavar='monomer/solvent', type=str, action='store', help='monomer or solvent order parameter.')
parser.add_argument('--png-name', dest='pn', metavar='imagename', type=str, action='store', help='Name of image file to be generated.', default="order_param")
parser.add_argument('--show-plot', dest='sp', action ='store_true' , help='Flag to include to see plot.') 
args = parser.parse_args() 

divnorm = matplotlib.colors.SymLogNorm (0.005, vmin=-0.1, vmax=0.1)

if __name__=="__main__":

    start = time.time() 
    #######################################
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes  () 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params (axis='x', labelsize=16)
    ax.tick_params (axis='y', labelsize=16)

    U_list    = aux.dir2U ( os.listdir (".") ) 
    # U_list    = ["U1"]
    PLOT_DICT = {} 

    starting_index = args.s
    ortn_file      = args.of
    dop            = args.dop
    show_plot_bool = args.sp
    i    = 0
    Tmax = 0 

    # instantiating pool 
    pool1 = multiprocessing.Pool ( processes=50 ) 
    pool2 = multiprocessing.Pool ( processes=5  )

    pool_list = [ pool1, pool2 ]

    for U in U_list:
        print ("Inside U = "+str(U)+"...", flush=True)
        temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(dop) ) )
        # temperatures = [0.01]
        # get the num_list for each temperature 
        master_temp_list = [] 
        master_num_list  = [] 
        ord_par1_dict     = {} 
        ord_par2_dict     = {} 
        ntraj_dict       = {} 

        # define stores 
        ord_par1_mean     = [] 

        for T in temperatures:
            num_list         = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend  ( num_list )
            master_temp_list.extend ( [T]* len(num_list) )
            ntraj_dict[T]    = len (num_list)
            ord_par1_dict [T] = []

        # start multiprocessing... keeping in mind that each node only has 96 cores
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list [0:50]
        mtemp_list_p2 = master_temp_list [50:100]
        mtemp_list_p3 = master_temp_list [100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3]
        ######################################################################
            
        mnum_list_p1  = master_num_list [0:50   ]
        mnum_list_p2  = master_num_list [50:100 ]
        mnum_list_p3  = master_num_list [100:105]
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3]

        # define a shitty dict
        shitty_dict = {0:0, 1:0, 2:1}

        for uidx in range(3):

            results = pool_list [ shitty_dict[uidx] ].starmap ( aux.obtain_order_parameter, zip( itertools.repeat(U), itertools.repeat(dop), \
                mtemp_list[uidx], itertools.repeat(ortn_file), mnum_list[uidx], itertools.repeat(starting_index) ) ) 
            print (results)
            print ("Pool has been closed. This pool has {} threads.".format ( len(results ) ), flush=True )

            for k in range ( len (mtemp_list[uidx]) ):
                # ord_par2_dict[ mtemp_list[uidx][k] ].append ( results[k][0] )
                # ord_par1_dict[ mtemp_list[uidx][k] ].append ( results[k][1] )
                ord_par1_dict[ mtemp_list[uidx][k] ].append ( results[k] )

            for T in np.unique ( mtemp_list [uidx] ):
                # ord_par2_mean.append ( np.mean ( ord_par2_dict[T] ) )
                # ord_par1_mean.append  ( np.mean  ( ord_par1_dict[T] ) )
                ord_par1_mean.append  ( np.mean  ( ord_par1_dict[T] ) )

        print (ord_par1_mean)
        PLOT_DICT [U] = np.asarray (ord_par1_mean)
    pool1.close()
    pool1.join () 

    pool2.close() 
    pool2.join ()

    i = 0 
    # colors = [ cm.seismic(x) for x in np.linspace(0, 1, len(U_list) ] 
    for U in U_list: 
        # print ("i=", i, ", len(U_list)=",len(U_list) )
        # print (PLOT_DICT[U])
        chi_a = aux.get_chi_cosolvent ( str(U)+"/geom_and_esurf.txt")[0]
        rgba_color = cm.PiYG_r (divnorm(chi_a))
        
        ax.errorbar ( temperatures, np.asarray(PLOT_DICT[U]), yerr=0, fmt='none', linestyle='-', elinewidth=1, capsize=2, linewidth=1, color='k', label='_nolegend_')
        ax.plot ( temperatures, (PLOT_DICT[U]), marker='o', markeredgecolor='k', linestyle='-', linewidth=2, c=rgba_color, label='_nolegend_', markersize=10)
        i += 1

    ##############################################################

    my_cmap = cm.PiYG_r
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=-0.1, vmax=0.1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [-0.1,0.1] )
    cbar.set_ticklabels( [-0.1, 0.1] ) 
    cbar.ax.tick_params(labelsize=14)
    ax.set_xscale('log')
    ax.set_ylim(0, 1)
    plt.savefig   ( args.pn + ".png", dpi=1000)
    
    if show_plot_bool:
        plt.show() 

    ##################################
    stop = time.time() 
    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

