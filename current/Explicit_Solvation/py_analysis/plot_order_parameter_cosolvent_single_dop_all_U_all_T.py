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

divnorm = matplotlib.colors.SymLogNorm ( 0.005, vmin=-0.1, vmax=0.1 )

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 

    # instantiate plt figure 
    plt.figure( figsize=(8,6) )

    i=0
    Tmax = []
    mm_max = [] 
    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        op_list = []
        op_err  = []
        op_mean = []
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
           
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 

                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "time_step"], engine='python', skiprows=args.s) 
                op_list.extend ( (df["ms1_tot"].values - df["ms2_tot"].values)/(df["ms1_tot"]+df["ms2_tot"]) )
                
            op_err.append( np.std(op_list)/np.sqrt(40) )
            op_mean.append( np.mean(op_list) ) 
        
        chi_1 = aux.get_chi_cosolvent ( str(U) + "/geom_and_esurf.txt")[1]
        rgba_color = cm.PiYG (divnorm (chi_2) ) 
        plt.plot ( temperatures, op_mean, marker='o', linestyle='-', color=rgba_color, markeredgecolor='k', linewidth=2, label='_nolegend_')
        # plt.errorbar(temperatures, np.asarray(op_mean), yerr=np.asarray(op_err), linestyle='-', elinewidth=1, capsize=0, color=cm.seismic(i/len(U_list)), fmt='o', markeredgecolor='k', label='_nolegend_')   
        i += 1
        # mm_max.append( np.max(mm_list) ) 
        print("done!", flush=True)
    
    contacts =  np.ones ( len(temperatures) ) 

    
    # # # # 
    my_cmap = cm.PiYG
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=-0.1, vmax=0.1) ) 
    ax = plt.axes() 
    

    plt.ylabel("$\\mu$", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    plt.xticks( temperatures, fontsize=16 )
    plt.yticks( fontsize=16 )
    ytick_list = np.unique( np.linspace (0, np.max(mm_max)+1, 20, dtype=int) )

    ytick_list = np.arange(0, ytick_list[-1], np.ceil(int(ytick_list[-1])/10 ) )  

    if not ( (np.max(mm_max)+1) in ytick_list):
        ytick_list = np.hstack ( (ytick_list, np.max(mm_max)+1 ) ) 

    # plt.yticks( ytick_list ) 

    cbar = plt.colorbar(sm, orientation='vertical')
    cbar.set_ticks( [-0.1, 0.1] )
    cbar.set_ticklabels( [-0.1, 0.1] )
    cbar.ax.set_ylabel ( "$ \chi ^1$ ", fontsize=16, rotation=270 ) 
    ax.set_xscale('log')
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    plt.yticks ( np.linspace(0,1,11) )
    # ax.yaxis.get_major_locator().set_params(integer=True)
    plt.savefig("DOP_"+str(args.dop)+"_cosolvent_order_param.png", dpi=1000)

    if args.sp:
        plt.show()
