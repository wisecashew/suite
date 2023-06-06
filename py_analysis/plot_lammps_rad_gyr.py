#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description = "Read a LAMMPS .avg file for gyration.")
parser.add_argument('-i', metavar='rad_gyr.avg', dest='i', action='store', help='enter address of xvg file', type=str)
args = parser.parse_args()

# get the original rdf

with open(args.i, 'r') as f:
    for line in f:
        if line.startswith('# TimeStep c_rog[0] c_rog[1] c_rog[2] c_rog[3]'): 
            break
    df = pd.read_csv(f, names=["t", "Rg", "Rgx", "Rgy", "Rgz"], delim_whitespace=True)
    f.close()

plt.plot(df["t"].values, df["Rg"].values) 
plt.xlabel("Time")
plt.ylabel("Radius of gyration") 
print( "Mean radius of gyration is {:0.2f}".format(np.mean(df["Rg"].values)) ) 

plt.show()
