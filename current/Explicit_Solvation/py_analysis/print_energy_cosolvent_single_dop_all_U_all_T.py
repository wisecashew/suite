#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('-T', dest='T', action='store', type=float, help='Provide a temperature.', default=100)
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide an energy surface.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 

args = parser.parse_args()

# divnorm = matplotlib.colors.SymLogNorm ( 0.005, vmin=-0.1, vmax=0.1 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    
    PLOT_DICT = {}

    i=0
    energy  = []
    Tlist   = []
    numlist = []

    num_list = np.unique ( aux.dir2nsim ( os.listdir ("." ) ) )

    for num in num_list: 
        df = pd.read_csv("energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
        energy.append (np.mean (df["energy"].values[-2000:]))

    print ( "Energy average = ", np.mean(energy) )

