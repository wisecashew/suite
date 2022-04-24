#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
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

i=0
for U in U_list:
    ms_list = np.array([])
    ms_mean = np.array([])
    ms_err  = np.array([]) 
    for temp in temperatures: 
        if U != "Uexcl":
            for idx in ["1", "2", "3", "4", "5"]:
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump"+idx+".txt", sep=' \| ', names=["energy", "mm", "ac", "nc", "time_step"], engine='python')
                ms_list = np.hstack ( (ms_list, np.mean( df["ac"].values[args.s:] + df["nc"].values[args.s:] ) ) )
            ms_mean = np.hstack  ( ( ms_mean, np.mean (ms_list) ) ) 
            ms_err  = np.hstack  ( ( ms_err , np.std  (ms_list)  / np.sqrt( 5 ) ) )  
        else:
            df = pd.read_csv( str(U) + "/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump.txt", sep=' \| ', names=["energy", "mm", "ac", "nc", "ts"], engine='python') 
            ms_mean = np.hstack ( ( ms_mean, np.mean( df["ac"].values[args.s:] + df["nc"].values[args.s:] ) ) )
            ms_err  = np.hstack ( ( ms_err , np.std ( df["ac"].values[args.s:] + df["nc"].values[args.s:] ) ) ) 

    if U == "Uexcl":
        plt.errorbar(temperatures, ms_mean, yerr=ms_err, linestyle='-', elinewidth=0.5, capsize=2, linewidth=3, alpha=0.5)  
    else:
        plt.errorbar(temperatures, ms_mean, yerr=ms_err, linestyle='-', elinewidth=1, linewidth=1.5, capsize=2, color=cm.copper(i/9))   
    
    i += 1

plt.legend(U_list, fontsize=12, loc='best')
plt.ylabel("Monomer-solvent contacts", fontsize=18)
plt.xlabel("Temperature", fontsize=18)
plt.xticks(np.arange(0,12,1))

plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.savefig("DOP_"+str(args.dop)+"_mscorr.png", dpi=1200)
plt.show()

