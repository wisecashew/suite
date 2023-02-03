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
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/py_analysis')
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
    lambda_T = np.array([-3, -1.4])
    df_T = pd.read_csv("TARGET/energydump_1.mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    avg_Nmm_T = np.mean  (df_T["mm_tot" ].values[-2000:])
    avg_Nms_T = np.mean  (df_T["ms1_tot"].values[-2000:]+df_T["ms2_tot"].values[-2000:])
    avg_T     = np.array ([avg_Nmm_T, avg_Nms_T])

    lambda_M = aux.get_energy_cg (args.m+"/geom_and_esurf.txt")
    
    df_M = pd.read_csv (args.m+"/energydump_1.mc",  sep=' \| ', names=["energy", "mm_tot", "ms_tot", "time_step"], engine='python', skiprows=0)
    avg_Nmm_M = np.mean (df_M["mm_tot"].values[-2000:])
    avg_Nms_M = np.mean (df_M["ms_tot"].values[-2000:])
    avg_M     = np.array([avg_Nmm_M, avg_Nms_M])
    avg_Nmm2_M = np.mean ((df_M["mm_tot"].values[-2000:])**2)
    avg_Nms2_M = np.mean ((df_M["ms_tot"].values[-2000:])**2)
    avg_M2     = np.array ([avg_Nmm2_M, avg_Nms2_M])
    
    correction = T*(avg_T - avg_M)/(avg_M2 - avg_M**2) 
    print (correction)
    lambda_M   = lambda_M - correction
    print ("New energies = ",lambda_M)
    f = open (args.m+"/delta_e.mc", 'w')
    f.write ("Emm = " + str(lambda_M[0]))
    f.write ("\n")
    f.write ("Ems = " + str(lambda_M[1]))
    f.write ("\n")
    f.write ("err = " + str( np.linalg.norm( correction ) ) )
    f.close()
    
