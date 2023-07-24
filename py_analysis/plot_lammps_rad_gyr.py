#!/home/satyend/.conda/envs/phase/bin/python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description = "Read a LAMMPS .avg file for gyration.")
parser.add_argument('-i', metavar='rad_gyr.avg', dest='i', action='store', help='enter address of xvg file', type=str)
args = parser.parse_args()

# get the original rdf

print ("Reading file.", flush=True)
with open(args.i, 'r') as f:
    for line in f:
        if line.startswith('# TimeStep c_rog[0] c_rog[1] c_rog[2] c_rog[3]'): 
            break
    df = pd.read_csv(f, names=["t", "Rg", "Rgx", "Rgy", "Rgz"], delim_whitespace=True)
    f.close()

print ("Done reading file.", flush=True)

print ("Plotting", flush=True)
rg_vals = df["Rg"].values[-10000:]
plt.hist(rg_vals, density=True, bins=50)
plt.savefig ("rg_dist", dpi=1200)

print( "Mean radius of gyration is {:0.2f}".format(np.mean(rg_vals)) ) 


