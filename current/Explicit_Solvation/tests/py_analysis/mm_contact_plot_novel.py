#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide potential energy surface.')
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl_vol', dest='ev', action='store', type=int, help='flag to include excluded volume forcefield.', default=0) 

args = parser.parse_args()  

temperatures = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0] 
U_list = [args.U]

plt.figure( figsize=(8,6) )

if args.ev==1:
    U_list.append("Uexcl")

i=0
for U in U_list:
    energy_list = []
    mm_list = []
    mm_err = []  
    for temp in temperatures: 
        df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/energydump.txt", sep=' \| ', names=["energy", "m_tot", "maa", "man", "mna", "mnn", "s_tot", "saa", "san", "sna", "snn", "ts"], engine='python', skiprows=1)
        mm_list.append(np.mean( df["m_tot"].values[args.s:] ) - (args.dop-1) )
        mm_err.append( np.std( df["m_tot"].values[args.s:] )  / np.sqrt(np.size( df["m_tot"].values[args.s::100] )) ) 

    if U == "Uexcl":
        plt.errorbar(temperatures, np.asarray(mm_list), yerr=np.asarray(mm_err), elinewidth=1, capsize=0, linewidth=3) # , linewidth=3) 
    else:
        plt.errorbar(temperatures, np.asarray(mm_list), yerr=np.asarray(mm_err), linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9))   
    i += 1
# datafile.close() 
print(mm_list)
plt.legend(U_list, fontsize=12, loc='best')
plt.ylabel("Average monomer-monomer contacts", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.xticks(np.arange(0,temperatures[-1]+2,1))

plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.savefig("DOP_"+str(args.dop)+"_mmcorr.png", dpi=1200)
plt.show()

