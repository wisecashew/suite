#!/usr/bin/env python3

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux
import os 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)

args = parser.parse_args()   

if __name__=="__main__":
    
    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir (".") )

    #instantiate plt figure
    plt.figure( figsize=(8,6) )

    # if we want the Uexcl contacts

    # instantiate some pertinent variables
    i=0
    Tmax = []
    ms_max = [] 
    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True )
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])

        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
            num_list = np.unique ( aux.dir2nsim (os.listdir ( str(U) + "/DOP_"+str(args.dop) + "/" + str(temp) ) ) )
            
            for num in num_list:
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)

                ms_list = np.hstack ( ( ms_list, ( df["ms_tot"].values ) - (args.dop-1) ) )
        
            ms_err  = np.hstack ( ( ms_err , np.std  ( ms_list ) / np.sqrt(40) ) ) 
            ms_mean = np.hstack ( ( ms_mean, np.mean ( ms_list ) ) )

        plt.errorbar(temperatures, np.asarray(ms_mean), yerr=np.asarray(ms_err), fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9), label='_nolegend_')   
        i += 1
        ms_max.append( np.max(ms_list) ) 
        print ("done!", flush=True)
    
    contacts = np.ones(len(temperatures))
    if args.ev:
        df = pd.read_csv ( "Uexcl/DOP_"+str(args.dop)+"/0.1/"+args.e, sep= ' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python')
        contacts = np.mean ( df["ms_tot"].values - (args.dop-1) ) * contacts
        plt.errorbar ( temperatures, contacts, yerr=np.std(df["ms_tot"].values)/np.sqrt(40), fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0 )
        plt.legend(["Athermal solvent"])

    my_cmap = cm.copper
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1))
    
    ax = plt.axes()

    plt.ylabel("Solvent contacts (SASA)", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    plt.xticks(temperatures, fontsize=12)

    cbar = plt.colorbar(sm, orientation='vertical')
    cbar.set_ticks( [0, 1] )
    cbar.set_ticklabels( ["Poorest", "Best"] )
    cbar.ax.set_ylabel ( "Quality of solvent", fontsize=14, rotation=270 ) 
    ax.set_xscale('log')
    plt.yticks(fontsize=12)
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    ax.yaxis.get_major_locator().set_params(integer=True)
    plt.savefig("DOP_"+str(args.dop)+"_mscorr.png", dpi=800)

    if args.sp:
        plt.show()

