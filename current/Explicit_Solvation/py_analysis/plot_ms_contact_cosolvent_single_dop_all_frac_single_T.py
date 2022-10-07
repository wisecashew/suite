#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.')
parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    
    # instantiate plt figure 
    plt.figure( figsize=(8,6) )
    
    PLOT_DICT  = {}
    ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    starting_index = args.s
    
    print ("T = ", args.T)
    temperatures = [float(elem) for elem in args.T]
    temperatures.sort() 

    ms_max = 25*2+(args.dop-2)*24
    for T in temperatures:
        
        frac_list    = [] 
        print ("Currently plotting out stuff in T = " + str(T) + "...", end=' ', flush=True)
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])

        for U in U_list: 
            frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
            ms_list = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=args.s) 
                ms_list = np.hstack ( (ms_list, df["ms1_tot"].values + df["ms2_tot"].values ) )
                
            ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(10) ) ) )
            ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )

        PLOT_DICT [T] = ms_mean
        ERROR_DICT[T] = ms_err
        # mm_max.append( np.max(ms_list) )
        print("done!", flush=True)
    
    i=0
    if len(temperatures) == 1:
        colors = cm.coolwarm (np.linspace(0,1,3)) 
        i = 1
    else:
        colors = cm.coolwarm (np.linspace(0,1,len(temperatures)))

    for T in temperatures:
        plt.errorbar ( frac_list, PLOT_DICT[T]/ ms_max, yerr=ERROR_DICT[T]/ms_max, linewidth=1, fmt='none', ecolor=colors[i], capsize=2, color='k', label="_nolabel_")
        plt.plot     ( frac_list, PLOT_DICT[T]/ ms_max, linestyle='-', marker='o',  markeredgecolor='k', linewidth=2, color=colors[i], label="_nolabel_", markersize=10)
        i += 1

    ####
    
    # my_cmap = cm.PiYG
    # sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0) ) 
    ax = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    
    # cbar = plt.colorbar(sm, orientation='vertical')
    # cbar.set_ticks( [0, 1.0] )
    # cbar.ax.tick_params (labelsize=14)
    # cbar.set_ticklabels( [0, 1.0] )
    # cbar.ax.set_ylabel ( "Strength of better solvent \n", fontsize=16, rotation=270 ) 
    # ax.set_xscale('log')
    # ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    # ax.yaxis.get_major_locator().set_params(integer=True)
    ax.set_xticks (np.linspace (0, 1, 11))
    ax.set_yticks (np.linspace (0.5,1,6))
    # plt.legend(U_list)
    plt.savefig   (args.pn + ".png", dpi=1000)

    if args.sp:
        plt.show()


