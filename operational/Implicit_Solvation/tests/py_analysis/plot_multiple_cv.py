#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting point for plotting.', default=100)
parser.add_argument('--excl_vol', dest='ev', action='store', type=int, help='flag to include excluded volume forcefield.', default=0)  

args = parser.parse_args()  


temperatures = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] 
U_list = ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9"]

if (args.ev==1):
    U_list.append("Uexcl") 
else:
    pass

dop = args.dop

plt.figure(1)
plt.title("Degree of polymerization = "+str(dop) ) 

for U in U_list:
    energy_list = []
    cv_list = []  
    for temp in temperatures: 
        df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/energydump.txt", sep=' \| ', names=["energy", "mm", "ac", "nc", "time_step"], engine='python')
        cv_list.append((np.mean( df["energy"].values[args.s:] **2 ) - (np.mean ( df["energy"].values[args.s:] ) **2 ))/(temp**2))
        energy_list.append(np.mean( df["energy"].values[args.s:] ) ) 
    plt.plot(temperatures, np.asarray(cv_list), linestyle='-.', marker='^', alpha=0.5)  

plt.legend(U_list)
plt.ylabel("heat capacity")
plt.xlabel("temperature")
plt.xticks(np.arange(0,11,1))
plt.show()

