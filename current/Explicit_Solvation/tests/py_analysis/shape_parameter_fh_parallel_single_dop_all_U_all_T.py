#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np 
import re 
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib.cm as cm 
import os
import argparse 
import aux
import multiprocessing
import itertools
import sys

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the average shape parameter for this trajectory.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number.', default=100)
parser.add_argument('--coords-file', dest='e', action='store', type=str, metavar='coords', help='Name of coords file to parse through.') 
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on a screen.', default=False)

args = parser.parse_args() 


def get_shape_parameter_from_one_traj ( U, DOP, T, num, coords_file, starting_index ):
    
    xlen = aux.edge_length (DOP) 
    ylen = aux.edge_length (DOP)
    zlen = aux.edge_length (DOP)

    filename = U + "/DOP_" + str(DOP) + "/" + str(T) + "/" + coords_file  + "_" + str(num) 
    master_dict = aux.get_pdict (filename, starting_index, DOP, xlen, ylen, zlen)
    
    shape_parameter = aux.get_shape_param ( master_dict, xlen, ylen, zlen ) 

    return shape_parameter 

##############################################################################################


if __name__=="__main__":
    
    U_list = aux.dir2U (os.listdir (".") )

    # instantiate plt figure 
    plt.figure ( figsize=(8,6) ) 
    ax  = plt.axes()
    dop = args.dop 
    xlen = aux.edge_length (dop)
    ylen = aux.edge_length (dop)
    zlen = aux.edge_length (dop) 

    # instantiate some pertinent variables 
    i = 0
    
    
    # instantiate pool...
    pool1 = multiprocessing.Pool ( processes=50 )
    pool2 = multiprocessing.Pool ( processes=5 ) 
    
    pool_list = [pool1, pool2] 

    ####################################################################
    i=0
    for U in U_list: 

        print ( "Currently plotting out stuff in U = " + str(U) + "...", flush=True ) 
        temperatures = aux.dir2float ( os.listdir ( U+"/DOP_"+str(dop) ) ) 
        shape_parameter = [] 
        
        shape_mean = [] 
        shape_std  = [] 
        # get num_list for each temperature 
        master_temp_list = []  # this is a temp list which goes [0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, ...] 
        master_num_list  = []  # this is a num list which goes [1,2,3,4,5, 1,2,3,4,5, ... ]
        shape_dict = {} # given a temperature, what is average shape factor 
        ntraj_dict = {} # given a temperature, how many trajectories exist 

        for T in temperatures: 
            num_list = aux.dir2nsim ( os.listdir ( U+"/DOP_"+str(dop)+"/"+str(T) ) ) 
            master_num_list.extend ( num_list ) 
            master_temp_list.extend ( [T]* len( num_list ) ) 
            ntraj_dict[T] = len (num_list) 
            shape_dict[T] = [] 
            
        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105] 
        mtemp_list = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1 = master_num_list[0:50]
        mnum_list_p2 = master_num_list[50:100]
        mnum_list_p3 = master_num_list[100:105]
        mnum_list = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 
        
        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1}

        for uidx in range(3):
            results = pool_list[ shitty_dict[uidx] ].starmap( get_shape_parameter_from_one_traj, \
                   zip(itertools.repeat(U), itertools.repeat(dop) , mtemp_list[uidx], \
                   mnum_list[uidx], \
                   itertools.repeat(args.e), itertools.repeat(args.s) ) ) 

            print ("Pool has done its job. This pool has {} threads.".format(len(results)), flush=True)

            for k in range ( len (mtemp_list[uidx] ) ):
                shape_dict[mtemp_list[uidx][k]].append ( results[k] ) 

            for T in np.unique ( mtemp_list[uidx] ): 
                shape_mean.append ( np.mean (shape_dict[T]) )
                shape_std.append  ( np.std  ( shape_dict[T] )/np.sqrt ( ntraj_dict[T] ) ) 

        ax.errorbar ( temperatures, shape_mean, yerr=shape_std, fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1 
       
        print ("done!")

    pool1.close() 
    pool1.join()

    pool2.close() 
    pool2.join()



    ################################################################
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )


    plt.ylabel ( "Shape parameter ($\\delta ^*$)", fontsize=18)
    plt.xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_xscale ("log")
    # plt.xticks ( temperatures) 
    cbar = plt.colorbar (sm, orientation='vertical') 
    cbar.set_ticks ( [0,1] ) 
    cbar.set_ticklabels ( ["Poorest", "Best"] ) 
    cbar.ax.set_ylabel ( "Quality of solvent", fontsize=16, rotation=270 )
    plt.yticks (np.arange(0,1,0.1), fontsize=12)
    plt.savefig ("DOP_"+str(dop)+"_shape_parameter.png", dpi=800)

    if args.sp:
        plt.show()

