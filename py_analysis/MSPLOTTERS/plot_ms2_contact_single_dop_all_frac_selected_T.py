#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
from matplotlib.ticker import StrMethodFormatter
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as tck
import argparse
import aux
import os


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.LogNorm ( vmin=0.01, vmax=100.0 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    fig = plt.figure(figsize=(4/2,3/2), constrained_layout=True)
    ax  = plt.axes()
    edge = aux.edge_length(args.dop)
    plt.rcParams["axes.labelweight"] = "bold"
    PLOT_DICT  = {}
    ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    
    print ("T = ", args.T)
    temperatures = [float(elem) for elem in args.T]
    temperatures.sort() 

    ms_max = 25*2+(args.dop-2)*24
    for T in temperatures:
        frac_list    = []
        npart_list   = []
        print ("Currently plotting out stuff in T = " + str(T) + "...", end=' ', flush=True)
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])

        for U in U_list: 
            frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
            nsol2 = int(np.floor(((edge**3)-dop)*frac_list[-1]))
            nsol1 = edge**3 - dop - nsol2
            npart_list.append ( (nsol1, nsol2) )
            mm_list = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
                mm_list = np.hstack ( (mm_list, np.mean(df["ms2_tot"].values[-2000:]) + np.mean(df["ms1_tot"].values[-2000:]) ) )
                
            mm_err  = np.hstack ( (mm_err , np.std (mm_list)/np.sqrt(len(num_list)) ) )
            mm_mean = np.hstack ( (mm_mean, np.mean(mm_list) ) )

        PLOT_DICT [T] = mm_mean/ms_max
        ERROR_DICT[T] = mm_err/ms_max

        print("done!", flush=True)
    
    i=0
    print (frac_list, flush=True)
    for T in temperatures:
        plt.errorbar ( frac_list, PLOT_DICT[T], yerr=ERROR_DICT[T], linewidth=1, fmt='none', ecolor='k', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( frac_list, PLOT_DICT[T], linestyle='-',  markeredgecolor='k', linewidth=3/2, color=cm.coolwarm(divnorm(T)), label="_nolabel_", markersize=4, marker='o')
        i += 1

    ###
    # write data points out for text later
    ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=8, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=8, direction="in", left="off", labelleft="on" )
    # ax.axhline (y=0, c='k', linewidth=1.0)
    # ax.yaxis.get_major_locator().set_params(integer=True)
    ax.minorticks_on()
    ax.set_xlim ((-0.05, 1.05))
    ax.set_ylim ((-0.05, 1.05))
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_xticklabels ([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks (np.arange (0, 1.1, 0.2))
    ax.set_yticklabels (ax.get_yticks(), weight='bold') 
    ax.set_xticklabels (ax.get_xticks(), weight='bold') 
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    plt.savefig   (args.pn, bbox_inches='tight', dpi=1200)

