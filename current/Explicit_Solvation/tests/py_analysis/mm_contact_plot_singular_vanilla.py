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

parser = argparse.ArgumentParser(description="Get the contacts for a given degree of polymerization for every energy surface.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide potential energy surface.')
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='flag to include excluded volume forcefield.', default=0) 

args = parser.parse_args()  

U_list = [args.U]

fig = plt.figure( figsize=(8,6) )


i=0
Tmax = []
for U in U_list:
    mm_mean = []
    mm_list = []
    mm_err  = [] 
   
    temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
    Tmax.append ( np.max(temperatures) )
    for temp in temperatures: 
        num_list = np.unique( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
        
        for num in num_list:
            df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mmaligned", "mmnaligned", "ms_tot", "msaligned", "msnaligned", "ts"], engine='python', skiprows=1)
            mm_list.append( df["mm_tot"].values[args.s:] - (args.dop-1) )

        mm_err.append( np.std ( mm_list ) / np.sqrt(5) ) 
        mm_mean.append( np.mean( mm_list ) )

    plt.errorbar(temperatures, np.asarray(mm_mean)/( args.dop*(args.dop-1)/2 ), yerr=np.asarray(mm_err)/ ( args.dop*(args.dop-1)/2 ), fmt='o', markeredgecolor='k', linestyle='-', color='darkgreen', elinewidth=1, capsize=0, label='_nolegend_' )   
    i += 1
    # print(mm_list)
# datafile.close() 
# print(mm_list)

if args.ev:
    mm_mean = [] 
    mm_err  = [] 
    mm_list = [] 

    temperatures = aux.dir2float ( os.listdir ( "Uexcl" + "/DOP_"+str(args.dop) ) ) 
    Tmax.append ( np.max(temperatures) ) 

    for temp in temperatures:
        num_list = np.unique ( aux.dir2nsim ( os.listdir ( "Uexcl/DOP_" + str(args.dop) + "/" + str(temp) ) ) ) 

        for num in num_list:
            df = pd.read_csv ( "Uexcl/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep='\|', names=["energy","mm_tot", "mmaligned", "mmnaligned", "ms_tot", "msaligned", "msnaligned", "ts"], engine='python', skiprows=1) 
            mm_list.append ( df["mm_tot"].values - (args.dop-1) ) 

        mm_err.append  ( np.std(mm_list)/np.sqrt(5) ) 
        mm_mean.append ( np.mean (mm_list) )

    plt.errorbar ( temperatures, np.asarray (mm_mean)/ (args.dop*(args.dop-1)/2), yerr=np.asarray(mm_err)/(args.dop*(args.dop-1)/2), fmt='^', markeredgecolor='k', linestyle='-.', elinewidth=1, capsize=0)
    

plt.legend(["Athermal solvent"], fontsize=12, loc='best')
plt.ylabel("$\\frac{\\langle M_C \\rangle}{(N\cdot (N-1)/2)}$", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.xticks(np.arange(0,np.max(Tmax)+2,1))
# ytick_list = np.unique( np.linspace (0, np.max(mm_mean)+1, 20, dtype=int) )
# plt.yticks( np.arange(0, ytick_list[-1], int(ytick_list[1]-ytick_list[0])+1 ) )
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.xscale('log')
# matplotlib.axis.Axis.set_major_formatter(formatter=matplotlib.ticker.ScalarFormatter())

plt.savefig("DOP_"+U+"_"+str(args.dop)+"_mmcorr.png", dpi=1200)
plt.show()

