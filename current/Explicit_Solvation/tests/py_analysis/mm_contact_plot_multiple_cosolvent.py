#!/usr/licensed/anaconda3/2020.7/bin/python

import re
import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux
import os 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl_vol', dest='ev', action='store', type=int, help='flag to include excluded volume forcefield.', default=0) 

args = parser.parse_args()  

U_list = aux.dir2U (os.listdir ('.') )# ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9"]
print("Ulist is: ", end=' ')
print(U_list)

plt.figure( figsize=(8,6) )

if args.ev==1:
    U_list.append("Uexcl")

i=0
Tmax = []
mm_max = [] 
for U in U_list:
    energy_list = []
    mm_list = []
    mm_err = []  
    temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
    Tmax.append ( np.max(temperatures) )
    for temp in temperatures: 
        # print(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump.txt")
        df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump.txt", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "time_step"], engine='python')
        
        mm_list.append(np.mean( df["mm_tot"].values[args.s:] ) - (args.dop-1) )
        mm_err.append( np.std( df["mm_tot"].values[args.s:] )  / np.sqrt(np.size( df["mm_tot"].values[args.s::100] )) ) 
        energy_list.append(np.mean( df["energy"].values[args.s:] ) ) 

    if U == "Uexcl":
        plt.errorbar(temperatures, np.asarray(mm_list), yerr=np.asarray(mm_err), elinewidth=1, capsize=0, linewidth=3)
    else:
        plt.errorbar(temperatures, np.asarray(mm_list), yerr=np.asarray(mm_err), linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9), marker='^')   
    i += 1
    mm_max.append( np.max(mm_list) ) 


plt.legend(["x=0", "x=0.1","x=0.25","x=0.36", "x=0.5", "x=0.6", "x=0.72", "x=0.8", "x=0.9","x=1.0"], fontsize=12, loc='best')
plt.ylabel("Average monomer-monomer contacts", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.xticks(np.arange(0,np.max(Tmax)+1,1))
plt.yticks(np.arange(0,np.max(mm_max)+1,1)) 

plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.savefig("DOP_"+str(args.dop)+"_mmcorr.png", dpi=1200)
plt.show()

