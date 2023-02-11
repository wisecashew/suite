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

def metro (x):
    # print (x)
    result = np.exp (-x)
    result[result>1]=1
    return result 

if __name__=="__main__":
    k = 1
    T = 1
    beta = 1/(k*T)
    energies_cg = aux.get_energy_cg (args.m+"/geom_and_esurf.txt")
    Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n  = aux.get_energy    ("TARGET/geom_and_esurf.txt")
    energies_t  = np.array ([Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n])
    df_T = pd.read_csv("TARGET/energydump_1.mc", sep='\|', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    Nmm_T         = df_T["mm_tot" ].values[-2000:]
    Nms_T         = df_T["ms1_tot"].values[-2000:] + df_T["ms2_tot"].values[-2000:]
    net_energy_T  = df_T["energy" ].values[-2000:] 
    net_energy_T  = net_energy_T - np.mean( net_energy_T )
    C_offset_T    = np.mean (net_energy_T)

    #######################################
    Emm_a_m, Emm_n_m, Ems1_a_m, Ems1_n_m, Ems2_a_m, Ems2_n_m, Es1s2_a_m, Es1s2_n_m  = aux.get_energy    (args.m+"/geom_and_esurf.txt")
    df_M = pd.read_csv (args.m+"/energydump_1.mc", sep='\|', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot", "ms1s2_aligned", "ms1s2_naligned", "timestep"], engine='python', skiprows=0)
    Nmm_M        = df_M["mm_tot"].values [-2000:]
    # print (Nmm_M)
    Nms_M        = df_M["ms1_tot"].values[-2000:] + df_M["ms2_tot"].values[-2000:]
    Nms1_a_M     = df_M["ms1_aligned" ].values[-2000:]
    Nms1_n_M     = df_M["ms1_naligned"].values[-2000:]
    net_energy_M = df_M["energy"].values [-2000:]
    net_energy_M = net_energy_M - np.mean( net_energy_M )
    C_offset_M   = np.mean (net_energy_M)

    U_CG   = Nmm_T*Emm_a_m + Nms_T*Ems1_a_m
    U_T    = net_energy_T
    dU     = np.mean ( (U_CG - np.mean(U_CG)) - (U_T - np.mean(U_T)) )

    ############ BAR routine 
    # numerator: sample target configuration
    beta = 1
    energy_M_in_T = Emm_a_m * Nmm_T + Ems1_a_m * Nms_T
    C_M           = np.mean (energy_M_in_T)
    energy_M_in_T = energy_M_in_T - np.mean (energy_M_in_T)
    energy_T_in_M = Emm_a * Nmm_M + Ems1_a_m * Ems1_a + Ems2_a_m * Ems1_n
    C_T           = np.mean (energy_T_in_M)
    C             = C_M - C_T
    energy_T_in_M = energy_T_in_M - np.mean (energy_T_in_M)
    numerator   = np.mean ( metro (beta*(energy_M_in_T - net_energy_T) ) )
    denominator = np.mean ( metro (beta*(energy_T_in_M - net_energy_M) ) )

    dfree_energy = -1/beta*np.log (numerator/denominator) + C
    print ("e_mm = " + str(Emm_a_m) + ", e_ms = " + str(Ems1_a_m))
    print ("free_energy = " + str(dfree_energy) )
    
    print ("<du/dlambda>_m - <du/dlambda>_t = ",np.mean (Nmm_M - Nmm_T))
    S_rel = dU - dfree_energy
    
    print ("S_rel = " + str(S_rel))
