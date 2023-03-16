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
parser.add_argument('--color-scheme', dest='cs', action='store', type=int, help='Type of coloring.')
parser.add_argument('--sim-style', dest='ss', action='store', type=str, help='flory or entropy.')
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()

if args.cs == 3:
    divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0.0, vmax=1.0 ) 
elif args.cs == 0:
    if args.ss == 'flory-huggins':
        divnorm = matplotlib.colors.SymLogNorm ( 0.0001, vmin=-0.2, vmax=0.1 ) # this is for flory-huggins  
    elif args.ss == 'entropy':
        divnorm = matplotlib.colors.SymLogNorm ( 0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 
    else:
        print ("Bad simstyle.", flush=True)
        exit()

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = ["U8", "U9", "U10"]
    U_list = ["U8", "U9", "U10"]
    plt.figure( figsize=(8,6) )
    
    PLOT_DICT = {}

    i=0
    Tmax = []
    ms_max = 25*2+(args.dop-2)*24
    # ms_max = -1
    for U in U_list:
        Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n = aux.get_energy (str(U)+"/geom_and_esurf.txt")
        print ("Emm_a = ",Emm_a,", Emm_n = ",Emm_n,", Ems1_a = ",Ems1_a,", Ems1_n = ",Ems1_n,", Ems2_a = ",Ems2_a,", Ems2_n = ",Ems2_n,", Es1s2_a = ",Es1s2_a,", Es1s2_n = ",Es1s2_n)
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        temperatures = [0.01,0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
        
        for temp in temperatures: 
            skip = 0
            ms_list = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
                f = df["mm_naligned"].values[-1500:]*Emm_n
                g = df["ms1_aligned"].values[-1500:]*Ems1_a
                ms_list = np.hstack ( (ms_list, np.mean( f*g ) - np.mean(f)*np.mean(g) ) )
                
            ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(len(num_list) ) ) ) )
            ms_mean = np.hstack ( (ms_mean, np.mean(ms_list) ) )

        PLOT_DICT[U] = ( ms_mean, ms_err ) 
        i += 1
        print("done!", flush=True)
    
    # get the maximum value of fluctuation
    for key in PLOT_DICT:
        if ms_max < np.max(PLOT_DICT[key][0]):
            ms_max = np.max(PLOT_DICT[key][0])

    print ("ms_max = ",ms_max)
    # ms_max = 1
    i=0
    f = open ("CV-COV-MMN-MSA-output.mc", 'w')
    ymin = 1
    chi_list = [-0.01, -0.1, -0.2]
    for U in U_list:
        rgba_color = cm.PiYG (divnorm (chi_list[i]) )
        plt.errorbar ( temperatures, PLOT_DICT[U][0]/(args.dop*np.array(temperatures)**2), yerr=PLOT_DICT[U][1]/np.array(temperatures)**2, linewidth=1, fmt='none', capsize=2, color='k', label="_nolabel_")
        plt.plot     ( temperatures, PLOT_DICT[U][0]/(args.dop*np.array(temperatures)**2), linestyle='-', marker='o',  markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
        if ymin > np.min( PLOT_DICT[U][0] ):
            ymin = np.min( PLOT_DICT[U][0] ) 
        i += 1
        f.write ("U = "+U+":\n")
        f.write ("T:  ")
        for t in temperatures:
            f.write ("{:>2.2f} ".format(t))
        f.write("\n")
        f.write ("Cv: ")
        for j in range(len(PLOT_DICT[U][0])):
            f.write ("{:>2.2e} ".format(PLOT_DICT[U][0][j]/temperatures[j]**2))
        f.write ("\n")

    f.close()


    ax = plt.axes ()
    ax.tick_params ( direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    
    ax.set_xscale('log')
    plt.gca().yaxis.set_major_formatter (StrMethodFormatter('{x:1.1e}'))
    ax.minorticks_on()
    ax.yaxis.set_minor_locator (matplotlib.ticker.AutoMinorLocator())
    plt.savefig (args.pn + ".png", bbox_inches='tight', dpi=1200)

