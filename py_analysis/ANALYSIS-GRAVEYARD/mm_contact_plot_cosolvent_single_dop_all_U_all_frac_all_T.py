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
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()  

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    PLOT_DICT = {} 
    # instantiate plt figure 
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params ( axis='x', labelsize=16 )
    ax.tick_params ( axis='y', labelsize=16 )
    i=0
    Tmax = []
    mm_max = -1
    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        mm_list = []
        mm_err  = []  
        mm_mean = []
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
           
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 

                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "time_step"], engine='python', skiprows=args.s) 
                mm_list.extend ( df["mm_tot"].values - (args.dop-1) )
                
            mm_err.append( np.std(mm_list)/np.sqrt(10000) )
            mm_mean.append( np.mean(mm_list) ) 

        # plt.errorbar(temperatures, np.asarray(mm_mean), yerr=np.asarray(mm_err), linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9), fmt='o', markeredgecolor='k', label='_nolegend_')   
        if mm_max < np.max (mm_mean):
            mm_max = np.max(mm_mean)

        PLOT_DICT[U] = (np.asarray(mm_mean), np.asarray(mm_err))
        print("done!", flush=True)
    
    i = 0
    for U in U_list:
        ax.errorbar ( temperatures, PLOT_DICT[U][0]/mm_max, yerr=PLOT_DICT[U][1]/mm_max, linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/len(U_list)), fmt='o', markeredgecolor='k', label='_nolegend_')
        i += 1 
    
    contacts =  np.ones ( len(temperatures) ) 

    if args.ev:
        df = pd.read_csv("Uexcl/DOP_"+str(args.dop)+"/0.1/energydump", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "time_step"], engine='python')
        contacts = np.mean ( df["mm_tot"].values - (args.dop-1) ) * contacts
        ax.errorbar ( temperatures, contacts/mm_max, yerr=0, fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0)
        # ax.legend (["Athermal solvent"], bbox_to_anchor=(90, 1), fontsize=12)
        ax.legend (["Athermal solvent"], loc='upper right', bbox_to_anchor=(1.1, 1.1), fontsize=12)
        # 

    
    # # # # 
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) ) 

    plt.ylabel("$\\langle M_C \\rangle / \\langle M_C \\rangle _{\\mathrm{max}}$", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    ytick_list = np.unique( np.linspace (0, np.max(mm_max)+1, 20, dtype=int) )

    ytick_list = np.arange(0, ytick_list[-1], np.ceil(int(ytick_list[-1])/10 ) )  

    if not ( (np.max(mm_max)+1) in ytick_list):
        ytick_list = np.hstack ( (ytick_list, np.max(mm_max)+1 ) ) 

    # plt.yticks( ytick_list ) 

    cbar = plt.colorbar(sm, orientation='vertical')
    cbar.set_ticks( [0, 1] )
    cbar.set_ticklabels( ["0", "0.2"] )
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ( "Fraction of cosolvent", fontsize=18, rotation=270 ) 
    ax.set_xscale('log')
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    ax.set_yticks ( np.linspace (0, 1, 11) )
    plt.savefig("DOP_"+str(args.dop)+"_cosolvent_multiple_mmcorr.png", dpi=1000)

    if args.sp:
        plt.show()


