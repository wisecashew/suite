#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import os 
import aux 

parser = argparse.ArgumentParser(description="Get the heat capacity for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting point for plotting.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0)  
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt') 
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.',default=False) 

args = parser.parse_args()  

U_list = aux.dir2U ( os.listdir(".") ) 

if (args.ev):
    U_list.append("Uexcl") 
else:
    pass

dop = args.dop

plt.figure(1)

i=1
for U in U_list:
    
    print ("Currently plotting out stuff in U = " + str(U) + "...", end =' ' )
    energy_list = np.asarray([])
    cv_list = np.asarray([])
    cv_err  = np.asarray([])
    cv_mean = np.asarray([])

    temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(args.dop) ) ) 

    for temp in temperatures: 
    
        num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(args.dop) + "/" + str(temp) ) ) )
        # print(num_list) 
        for num in num_list:
            df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot","ms_aligned", "ms_naligned", "time_step"], engine='python')
            energy_list = df["energy"].values[args.s:] 
            # print (energy_list)
            cv = (np.mean( (energy_list) ** 2 ) - np.mean ( energy_list )**2)/((temp**2)*args.dop) 
            cv_list = np.hstack ( (cv_list, cv ) ) 
        
        # print (cv_mean)
        # print (cv_err)
        cv_mean = np.hstack ( ( cv_mean, np.mean ( cv_list ) ) ) 
        cv_err  = np.hstack ( ( cv_err , np.std ( cv_list )/np.sqrt(5) ) ) 
        
    if U == "Uexcl":
        plt.errorbar ( temperatures, np.asarray ( cv_mean ), yerr = np.asarray(cv_err), fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=3) 
    else:
        plt.errorbar ( temperatures, np.asarray ( cv_mean ), yerr = np.asarray (cv_err), fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9), label='_nolegend_' ) 

    i += 1
    print ("done!") 


my_cmp = cm.copper 
sm = plt.cm.ScalarMappable (cmap=my_cmp, norm=plt.Normalize(vmin=0, vmax=1))

cbar = plt.colorbar ( sm, orientation='vertical' )
cbar.set_ticks ( [0, 1] ) 
cbar.set_ticklabels ( ["Poor solvent", "Good solvent"] ) 
cbar.ax.set_ylabel ( "Solvent quality", fontsize=18, rotation=270 ) 
plt.legend(["Athermal solvent"])
plt.ylabel("$C_v/N$", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.xticks(np.arange(0,11,1))
plt.savefig ( "DOP_"+str(args.dop)+"_CV.png", dpi=1200)

if args.sp:
    plt.show()



