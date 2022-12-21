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
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number (not index or line number in file).', default=100)
parser.add_argument('-T', metavar='T', dest='T', action='store', nargs='+', type=float, help='Enter temperatures you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('-ymax', metavar='ymax', type=float, dest='ymax', action='store', default=2.8, help='Request upper limit of plot.')
parser.add_argument('--coords', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--png-name', dest='pn', metavar='imagename', action='store', type=str, help='Name of image file', default='rg_plot')
args = parser.parse_args() 

divnorm = matplotlib.colors.LogNorm ( vmin=0.01, vmax=100 ) 


def get_starting_ind ( U, T, num, dop, dumpfile):
    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/" + dumpfile + "_" + str(num) + ".mc"
    df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    L = len(df["energy"])
    return int(df["time_step"].values[L-2000])


def infiltrate_coords_get_rg ( U, T, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num)+".mc" 
    edge = aux.edge_length (dop)
    master_dict = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
    rg = aux.get_Rg(master_dict, edge, edge, edge) 
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
    
    U_list = aux.dir2U ( os.listdir (".") )
    print (U_list)
    PLOT_DICT = {} 
    dop            = args.dop
    coords_files   = args.c
    starting_index = args.s
    
    print ("T = ", args.T)
    temperatures = [float(elem) for elem in args.T] 
    temperatures.sort() 
    ######
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    
    # rg_max = 1/(6**0.5) * (dop**0.5)
    rg_max = (1+np.sqrt(2)+np.sqrt(3))/(3*6**0.5) * (dop**0.57) 
    # instantiating pool
    nproc = args.nproc
    pool1 = multiprocessing.Pool ( processes=nproc )# len(num_list)) 
    
    pool_list = [pool1] # , pool2]
    
    f = open("RG_DATA_"+str(dop)+".mc", "w") 

    for T in temperatures:
        frac_list = [] 
        f.write ( "T = " + str(T) + ":\n" )
        print("Inside T = " + str(T) + ", and N = " + str(dop) + "...", flush=True )
        rg_mean = [] 
        rg_std  = [] 
        
        # get num_list for each temperature 
        master_U_list     = []
        master_num_list   = []
        master_index_list = []
        rg_dict    = {}
        ntraj_dict = {}
        for U in U_list:
            frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt") )
            num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            for num in num_list:
                master_index_list.append (get_starting_ind (U, T, num, dop, "energydump") )
            master_U_list.extend ( [U]*len( num_list ) )
            ntraj_dict[U] = len ( num_list )
            rg_dict[U] = []

        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        
        idx_range = len (master_num_list)//nproc + 1
        for u_idx in range (idx_range):
            if u_idx == idx_range-1:
                results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( master_U_list[u_idx*nproc:], itertools.repeat(T), master_num_list [u_idx*nproc:], itertools.repeat (dop), itertools.repeat (coords_files), master_index_list[u_idx*nproc:] ) )
            else:
                results = pool_list[ 0 ] .starmap ( infiltrate_coords_get_rg, zip( master_U_list[(u_idx)*nproc:(u_idx+1)*nproc], itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(coords_files), master_index_list[u_idx*nproc:(u_idx+1)*nproc] ) )

            print ("Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

            for k in range( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
                rg_dict[master_U_list[u_idx*nproc + k]].append( results[k] )
        
        sorted_U_list = list(np.unique (master_U_list))
        sorted_U_list.sort(key=mysorter_f)
        # print ("sU list: ", sorted_U_list)
        for U in sorted_U_list:
            # print ("U = ", U, ", rg_dict[", U, "] = ", rg_dict[U])
            rg_mean.append( np.mean ( rg_dict[U] ) )
            rg_std.append ( np.std  ( rg_dict[U] ) / np.sqrt(master_U_list.count(U) ) )
        
        # print (rg_mean)
        # if rg_max < np.max (rg_mean):
        #     rg_max = np.max(rg_mean)
        
        PLOT_DICT [T] = (np.asarray(rg_mean), np.asarray(rg_std))
         
        f.write("Rg^2: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("x: ") 
        for elem in frac_list: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()
    f.close() 
    print ("rg_max = ", rg_max)
    i=0
    if len(temperatures) == 1:
        colors = cm.coolwarm (np.linspace(0,1,3)) 
        i = 1
    else:
        colors = cm.coolwarm (np.linspace(0,1,len(temperatures)))
    
    for T in temperatures:
        print (temperatures)
        ax.errorbar ( frac_list, PLOT_DICT[T][0]/rg_max, yerr= PLOT_DICT[T][1]/rg_max, linewidth=1, capsize=2, color='k', fmt='none', label='_nolegend_')
        ax.plot     ( frac_list, PLOT_DICT[T][0]/rg_max, marker='o', markeredgecolor='k', linestyle='-', linewidth=3, c=cm.coolwarm(divnorm(T)), label='_nolegend_', markersize=10 ) 
        i += 1

    #########
    # plt.rcParams['text.usetex'] = True
    my_cmap = cm.coolwarm
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0) )
    ax.minorticks_on()
    ax.set_xlim((-0.05, 1.05))
    ax.set_ylim(0.2, 1.6)
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_yticks (np.arange(0.3, 1.5, 0.3))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.2f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    plt.savefig   ( args.pn+".png", bbox_inches='tight', dpi=1200)
    
    ##################################
    stop = time.time() 

    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

