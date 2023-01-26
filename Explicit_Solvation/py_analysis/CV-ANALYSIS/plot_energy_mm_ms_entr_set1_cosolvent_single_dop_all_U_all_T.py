#!/usr/licensed/anaconda3/2020.7/bin/python

import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis")
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
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    # U_list = aux.dir2U ( os.listdir(".") ) 
    U_list = ["U11"]
    plt.figure( figsize=(8,6) )

    PLOT_DICT = {}

    i=0
    Tmax = []
    ms_max = 25*2+(args.dop-2)*24
    for U in U_list:
        Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n = aux.get_energy (str(U)+"/geom_and_esurf.txt")
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])
        # temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        temperatures = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
        
        for temp in temperatures: 
            skip = 0
            ms_list = np.asarray ([])
            mm_list = np.asarray ([])
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
                f = df["mm_aligned" ].values[-2000:]*Emm_a  + df["mm_naligned" ].values[-2000:]*Emm_n
                g = df["ms1_aligned"].values[-2000:]*Ems1_a + df["ms1_naligned"].values[-2000:]*Ems1_n
                ms_list = np.hstack ( ( ms_list, np.mean ( g ) ) )
                mm_list = np.hstack ( ( mm_list, np.mean ( f ) ) )
                
            ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(30) ) ) )
            ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )
            mm_err  = np.hstack ( (mm_err,  (np.std(mm_list)/np.sqrt(30) ) ) )
            mm_mean = np.hstack ( (mm_mean, np.mean(mm_list) ) )

        PLOT_DICT[U] = ( ms_mean, ms_err, mm_mean, mm_err ) 
        i += 1
        # mm_max.append( np.max(ms_list) ) 
        print("done!", flush=True)
    
    i=0
    ymin = 1
    ms_max = 1
    chi_list = [0.1, 0.05, 0.01, 0.001, 0, -0.001, -0.01, -0.1, -0.2]
    for U in U_list:
        rgba_color = cm.PiYG (divnorm (chi_list[i]) )
        plt.errorbar ( temperatures, PLOT_DICT[U][0] , yerr=PLOT_DICT[U][1], linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( temperatures, PLOT_DICT[U][0] , linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color='steelblue', label="_nolabel_", markersize=10)
        plt.errorbar ( temperatures, PLOT_DICT[U][2] , yerr=PLOT_DICT[U][3], linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( temperatures, PLOT_DICT[U][2] , linestyle='-', marker='^',  markeredgecolor='k', linewidth=3, color='salmon', label="_nolabel_", markersize=10)
        i += 1

    # plot excluded volume
    ax = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    
    ax.set_xscale('log')
    # ax.set_ylim((0.0 , 1.06))
    plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.0f}'))
    # ax.set_yticks (np.linspace (0.0,1,6))
    ax.minorticks_on()
    # ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
    plt.savefig (args.pn + ".png", bbox_inches='tight', dpi=1200)


