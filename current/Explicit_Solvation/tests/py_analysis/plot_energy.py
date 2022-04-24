#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide name of energy surface.') 
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting point for plotting.', default=100)

args = parser.parse_args()  


temperatures = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] 
dop = args.dop 
U   = args.U   

plt.figure(1)
plt.title("$U$ = "+str(args.U) ) 


energy_list = []
energy_std  = []
for temp in temperatures: 
    df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/energydump3.txt", sep=' \| ', names=["energy", "mm", "ac", "nc", "time_step"], engine='python')
    # cv_list.append((np.mean( df["energy"].values[args.s:] **2 ) - (np.mean ( df["energy"].values[args.s:] ) **2 ))/(temp**2))
    energy_list.append(np.mean( df["energy"].values[args.s:] ) ) 
    energy_std.append (np.std ( df["energy"].values[args.s:] ) )

plt.errorbar(temperatures, np.asarray(energy_list), yerr=energy_std, linestyle='-.', marker='^')  

# plt.legend(Np_list)
plt.ylabel("heat capacity")
plt.xlabel("temperature")
plt.xticks(np.arange(0,11,0.5))
plt.show()

