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
    
    # instantiate plt figure 
    plt.figure( figsize=(8,6) )
    
    PLOT_DICT  = {}
    ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    # starting_index = args.s
    
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
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
                ms_list = np.hstack ( (ms_list, np.mean(df["ms1_tot"].values[-2000:] + df["ms2_tot"].values[-2000:]) ) )
                
            ms_err  = np.hstack ( (ms_err,  np.std(ms_list)/np.sqrt(30) ) ) 
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
        colors = cm.coolwarm (np.array(temperatures) )# cm.coolwarm ( np.logspace ( int ( np.log (temperatures[0] ) ),int ( np.log(temperatures[-1] ) ), len ( temperatures ) ) )

    for T in temperatures:
        plt.errorbar ( frac_list, PLOT_DICT[T]/ ms_max, yerr=ERROR_DICT[T]/ms_max, linewidth=1, fmt='none', ecolor='k', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( frac_list, PLOT_DICT[T]/ ms_max, linestyle='-',  markeredgecolor='k', linewidth=3, color=cm.coolwarm(divnorm(T)), label="_nolabel_", markersize=10, marker='o')
        i += 1

    ###
    ax = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.minorticks_on()
    ax.set_xlim((-0.05, 1.05))
    ax.set_ylim((0.48 , 1.02))
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks (np.linspace (0.50,1,6))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.2f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    

    plt.savefig   (args.pn + ".png", bbox_inches='tight', dpi=1200)


