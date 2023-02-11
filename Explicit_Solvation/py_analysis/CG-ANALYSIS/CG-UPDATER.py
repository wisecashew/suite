#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse 
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/lattice_md/Explicit_Solvation/py_analysis')
import aux 
import os 

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument ("--model", dest='m', type=str, help="Select model directory.")
args = parser.parse_args()

if __name__=="__main__":

    # get the entire list of potential energy surfaces. 
    # step 1: find <dU_M/dlambdai>_T. 
    # T is the target ensemble: the real ensemble. 
    # M is the coarse-grained potential energy function. 
    # two lambdas: E_mm, E_ms. 
    # step 1 involves finding <N_mm> and <N_ms>
    # step 1 (recast): find <N_mm> and <N_ms> 
    T = 1
    df_T = pd.read_csv("TARGET/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    avg_Nmm_T = np.mean  (df_T["mm_tot" ].values[-2000:])
    avg_T     = avg_Nmm_T 
    print ("target average mm contacts = " + str(avg_Nmm_T))
    print ("target average contacts = ",avg_T)
   
    energies = aux.get_energy_cg (args.m+"/geom_and_esurf.txt")
    lambda_M = energies[0] - 2*energies[1]
    print ("parameters to update = ",lambda_M)
    
    df_M = pd.read_csv (args.m+"/energydump_1.mc",  sep=' \| ', names=["energy", "mm_tot", "ms_tot", "time_step"], engine='python', skiprows=0)

    avg_Nmm_M = np.mean (df_M["mm_tot"].values[-2000:])
    avg_M     = avg_Nmm_M

    print ("model average mm contacts = " + str(avg_Nmm_M))
    print ("model average contacts = ",avg_M)
    print ("model2 average contacts = ",avg_M**2)

    avg_Nmm2_M = np.mean  ((df_M["mm_tot"].values[-2000:])**2)
    avg_M2     = avg_Nmm2_M

    print ("model average mm contacts2 = " + str(avg_Nmm2_M))
    print ("model average contacts2 = ",avg_M2)

    correction = T*(avg_T - avg_M)/(avg_M2 - avg_M**2)
    print ("correction from initial guess = ",correction)
    lambda_M   = lambda_M - correction
    print ("New energies = ",lambda_M)
    f = open (args.m+"/delta_e.mc", 'w')
    f.write ("Emm = " + str( energies[0] )) # lambda_M + 2*energies[1]  ))
    f.write ("\n")
    f.write ("Ems = " + str( (energies[0]-lambda_M)/2 ) )
    f.write ("\n")
    f.write ("err = " + str( np.linalg.norm( correction ) ) )
    f.close()
    
