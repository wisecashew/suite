#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
# matplotlib.use('Agg')
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

divnorm = matplotlib.colors.SymLogNorm(0.005, vmin=-0.1, vmax=0.1) 

if __name__=="__main__":
    
    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir (".") )
    #instantiate plt figure
    plt.figure( figsize=(8,6) )
    ax = plt.axes() 
    ax.tick_params ( axis='x', labelsize=12, direction="in", left="off", labelleft="on" )
    ax.tick_params ( axis='y', labelsize=12, direction="in", left="off", labelleft="on" )
    # instantiate some pertinent variables
    i=0
    Tmax = []
    # mm_max = [] 
    mm_max  = -1
    for U in U_list:
        if U == "U1":
            s = 50000
            
        else:
            s = args.s
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True )
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
            num_list = np.unique ( aux.dir2nsim (os.listdir ( str(U) + "/DOP_"+str(args.dop) + "/" + str(temp) ) ) )
            
            for num in num_list:
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=s)
                mm_list = np.hstack ( ( mm_list, ( df["mm_tot"].values ) - (args.dop-1) ) )

            if U == "U1" and temp == 0.1:
                mm_max = np.max( mm_list )
            mm_err  = np.hstack ( ( mm_err , np.std  ( mm_list ) / np.sqrt(40) ) ) 
            mm_mean = np.hstack ( ( mm_mean, np.mean ( mm_list ) ) )
        
        chi_a = aux.get_chi_entropy(str(U)+"/geom_and_esurf.txt")[0]
        rgba_color = cm.coolwarm(divnorm(chi_a))
        # print ("rgba color is ", rgba_color)
        ax.plot (temperatures, np.asarray(mm_mean)/mm_max, marker='o', markeredgecolor='k', linestyle='-', linewidth=2, c=rgba_color, label='_nolegend_')
        i += 1
        # mm_max.append( np.max(mm_list) )
        print ("done!", flush=True)
        # print ( mm_list )
    
    
    contacts = np.ones ( len (temperatures) )
    
    # my_cmap = cm.coolwarm
    # sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=-0.1, vmax=0.1))
    
    ax = plt.axes()

    plt.ylabel("$\\langle M_C \\rangle/\\langle M_C \\rangle _{\\mathrm{max}}$", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    # plt.xticks( temperatures, fontsize=12 )
    # cbar = plt.colorbar(sm, orientation='vertical')
    # cbar.set_ticks( [-0.1, 0.1] )
    #cbar.set_ticklabels( ["-0.1", "0.1"] )
    # cbar.ax.tick_params(labelsize=14)
    # cbar.ax.set_ylabel ( "$\chi ^a$", fontsize=18, rotation=270, labelpad=0 ) 
    ax.set_xscale('log')
    # ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    # ax.set_yticks ( np.linspace (0, 1, 11) )
    plt.savefig("DOP_"+str(args.dop)+"_multiple_mmcorr.png", dpi=1000)

    if args.sp:
        plt.show()

