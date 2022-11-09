#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--color-scheme', dest='cs', action='store', type=int, help='Type of coloring.')
parser.add_argument('--sim-style', dest='ss', action='store', type=str, help='flory or entropy.')
parser.add_argument('--path-to-excl', dest='pte', metavar='add/ress', action='store', type=str, help='Path to file.')
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

if args.cs == 3:
    divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0.0, vmax=1.0 ) 
elif args.cs == 0:
    # divnorm = matplotlib.colors.SymLogNorm ( 0.0001, vmin=-0.2, vmax=0.1 ) # this is for flory-huggins  
    # divnorm = matplotlib.colors.SymLogNorm ( 0.05, vmin=-0.2, vmax=0.1 ) # this is for entropy 
    divnorm = matplotlib.colors.SymLogNorm ( 0.05, vmin=-1000, vmax=1000 )

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    U_list.remove("U8")
    U_list.remove("U7")
    # U_list = ["U1"]
    plt.figure( figsize=(8,6) )
    
    PLOT_DICT = {}

    i=0
    Tmax = []
    # ms_max = 25*2+(args.dop-2)*24
    ms_max = -1
    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        temperatures = [0.01, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
        
        for temp in temperatures: 
            skip = args.s
            if temp < 2:
               skip = 7000
            elif temp == 2.5 or temp == 5.0:
                skip = 9000
            #else:
            #    skip = args.s
            ms_list = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
                f = df["ms1_aligned"].values/df["ms1_tot"].values
                ms_list = np.hstack ( (ms_list, np.mean( f**2 )-(np.mean(f) )**2  ) )
                
            ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(30) ) ) )
            ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )

        PLOT_DICT[U] = ( ms_mean, ms_err ) 
        i += 1
        # mm_max.append( np.max(ms_list) ) 
        print("done!", flush=True)
    
    # get the maximum value of fluctuation
    for key in PLOT_DICT:
        if ms_max < np.max(PLOT_DICT[key][0]):
            ms_max = np.max(PLOT_DICT[key][0])

    print ("ms_max = ",ms_max)
    i=0
    ms_max = 1
    ymin = 1
    for U in U_list:
        chi_1 = aux.get_chi_cosolvent ( str(U)+"/geom_and_esurf.txt" )[args.cs]
        print ("chi_1 = ", chi_1)
        if args.cs == 0: 
            rgba_color = cm.PiYG (divnorm (chi_1) )
        else:
            rgba_color = cm.PiYG(divnorm (chi_1))
        plt.errorbar ( temperatures, PLOT_DICT[U][0] / ms_max, yerr=PLOT_DICT[U][1]/ms_max, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( temperatures, PLOT_DICT[U][0] / ms_max, linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
        if ymin > np.min( PLOT_DICT[U][0]/ms_max ):
            ymin = np.min( PLOT_DICT[U][0]/ms_max ) 
        i += 1

    # plot excluded volume
    contacts =  np.ones ( len(temperatures) ) 
    if args.ev:
        print ("In U = Uexcl...", flush=True)
        df = pd.read_csv(args.pte+"/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "s1s2_tot", "s1s2_aligned", "s1s2_naligned", "time_step"], engine='python')
        contacts = np.mean ( df["ms1_tot"].values + df["ms2_tot"].values ) * contacts
        plt.errorbar ( temperatures, contacts/ms_max, yerr=0, fmt='-', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, markersize=10, linewidth=3)

    
    if args.cs == 0:
        my_cmap = cm.PiYG
        sm = plt.cm.ScalarMappable( cmap=my_cmap, norm=plt.Normalize(vmin=-0.2, vmax=0.1) )
    else:
        my_cmap = cm.PiYG
        sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=-0.5, vmax=0.5 ) )
    ax = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    
    cbar = plt.colorbar (sm, orientation='vertical', format='%.0e')
    cbar.set_ticks ( [-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1] )
    cbar.ax.tick_params (labelsize=14)
    cbar.set_ticklabels( ["-$10^{-1}$", "-$10^{-2}$","-$10^{-3}$",  0, "$10^{-3}$", "$10^{-2}$", "$10^{-1}$" ] )
    ax.set_xscale('log')
    # ax.set_ylim((0.48 , 1.02))
    plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.0e}'))
    # ax.set_yticks (np.linspace (0.50,1,6))
    ax.minorticks_on()
    ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
    plt.savefig (args.pn + ".png", dpi=1000)

    if args.sp:
        plt.show()


