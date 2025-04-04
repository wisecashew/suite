#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl_vol', dest='ev', action='store', type=int, help='flag to include excluded volume forcefield.', default=0)

args = parser.parse_args()  


temperatures = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] 
U_list =["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9"]

plt.figure()
plt.title("DOP = "+str(args.dop), fontsize=24 ) 

if args.ev==1:
    U_list.append("Uexcl") 


for U in U_list:
    energy_list = []
    ms_list = []
    ms_err = []  
    for temp in temperatures: 
        df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump.txt", sep=' \| ', names=["energy", "mm", "ac", "nc", "time_step"], engine='python')
        ms_list.append(np.mean( df["ac"].values[args.s:] + df["nc"].values[args.s:] ) - (args.dop-1) )
        ms_err.append( np.std( df["ac"].values[args.s:] + df["nc"].values[args.s:] )  / np.sqrt(np.size( df["mm"].values[args.s:] )) ) 
        energy_list.append(np.mean( df["energy"].values[args.s:] ) ) 
    if U == "Uexcl":
        plt.errorbar(temperatures, np.asarray(ms_list), yerr=np.asarray(ms_err), linestyle='-', elinewidth=1, capsize=2, linewidth=3)  
    else:
        plt.errorbar(temperatures, np.asarray(ms_list), yerr=np.asarray(ms_err), linestyle='-.', elinewidth=1, capsize=2)   


plt.legend(U_list, fontsize=12, loc='best')
plt.ylabel("Monomer-solvent contacts", fontsize=18)
plt.xlabel("Temperature", fontsize=18)
plt.xticks(np.arange(0,12,1))

plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.savefig("DOP_"+str(args.dop)+"_mscorr.png", dpi=1200)
# plt.show()

