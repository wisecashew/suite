#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import re 
import argparse 
import os 
import aux 


''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide potential energy surface.')
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0) 

args = parser.parse_args()  

U_list = [args.U]

plt.figure( figsize=(8,6) )

if args.ev:
    U_list.append("Uexcl")

i=0
Tmax = []
for U in U_list:
    mm_mean = np.asarray([])
    mm_list = np.asarray([])
    mm_err  = np.asarray([]) 
   
    temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
    Tmax.append ( np.max(temperatures) )
    for temp in temperatures: 
        num_list = np.unique( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
        
        for num in num_list:
            df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mmaligned", "mmnaligned", "ms_tot", "msaligned", "msnaligned", "ts"], engine='python', skiprows=1)
            mm_list = np.hstack( (mm_list, ( df["mm_tot"].values[args.s:] ) - (args.dop-1) ) )

        mm_err  = np.hstack( (mm_err , np.std ( mm_list ) / np.sqrt(5) ) ) 
        mm_mean = np.hstack( (mm_mean, np.mean( mm_list ) ) ) 

    if U == "Uexcl":
        plt.errorbar(temperatures, np.asarray(mm_mean), yerr=np.asarray(mm_err), elinewidth=1, capsize=0, linewidth=3) # , linewidth=3) 
    else:
        plt.errorbar(temperatures, np.asarray(mm_mean), yerr=np.asarray(mm_err), linestyle='-.', elinewidth=1, capsize=0, color=cm.copper(i/9), label='_nolegend_' )   
    i += 1
    # print(mm_list)
# datafile.close() 
# print(mm_list)
plt.legend(U_list, fontsize=12, loc='best')
plt.ylabel("Average monomer-monomer contacts", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.xticks(np.arange(0,np.max(Tmax)+2,1))
ytick_list = np.unique( np.linspace (0, np.max(mm_mean)+1, 20, dtype=int) )
plt.yticks( np.arange(0, ytick_list[-1], int(ytick_list[1]-ytick_list[0])+1 ) )
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.xscale('log')
plt.savefig("DOP_"+str(args.dop)+"_mmcorr.png", dpi=1200)
plt.show()

