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
parser.add_argument('--png-name', dest='pn', metavar='imagename', action='store', type=str, help='Name of image file.', default="mm_plot")
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.005, vmin=-0.1, vmax=0.1 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    
    # instantiate plt figure 
    plt.figure( figsize=(8,6) )
    
    PLOT_DICT = {}

    i=0
    Tmax = []
    mm_max = 28-(args.dop-1)
    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
            mm_list = np.asarray([])
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 
                if temp == 0.01:
                    skip = 1500
                else:
                    skip = args.s

                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned","s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], engine='python', skiprows=skip)
                mm_list = np.hstack( (mm_list, df["mm_tot"].values - (args.dop-1) ) )
                # print (df)

            mm_err  = np.hstack((mm_err, np.std(mm_list)/np.sqrt(80)))
            mm_mean = np.hstack((mm_mean, np.mean(mm_list)))

        PLOT_DICT[U] = ( mm_mean, mm_err )
        # print (mm_mean)
        print("done!", flush=True)
    
    i=0
    for U in U_list:
        chi_1 = aux.get_chi_cosolvent ( str(U)+"/geom_and_esurf.txt" )[0]
        rgba_color = cm.PiYG(divnorm (chi_1))
        plt.errorbar ( temperatures, PLOT_DICT[U][0]/mm_max, yerr=PLOT_DICT[U][1]/mm_max, linewidth=1, fmt='none', capsize=2, color='k', label='_nolegend_')
        plt.plot ( temperatures, PLOT_DICT[U][0]/mm_max, linestyle='-', marker='o', color=rgba_color,  markeredgecolor='k', linewidth=2, label='_nolegend_', markersize=10)
        i += 1

    # plot excluded volume
    contacts =  np.ones ( len(temperatures) ) 
    if args.ev:
        df = pd.read_csv("Uexcl/DOP_"+str(args.dop)+"/0.1/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned","s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], engine='python', skiprows=args.s)
        contacts = np.mean ( df["mm_tot"].values - (args.dop-1) ) * contacts
        # print (contacts)
        plt.errorbar ( temperatures, contacts/mm_max, yerr=0, fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0)

    
    # # # # 
    my_cmap = cm.PiYG
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=-0.1, vmax=0.1) ) 
    ax = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    
    cbar = plt.colorbar(sm, orientation='vertical')
    cbar.set_ticks( [-0.1, 0.1] )
    cbar.ax.tick_params (labelsize=14)
    cbar.set_ticklabels( [-0.1, 0.1] )
    ax.set_xscale('log')
    ax.set_yticks (np.linspace (0,1,11))
    plt.savefig(args.pn+".png", dpi=1000)

    if args.sp:
        plt.show()


