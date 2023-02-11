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

    #######################################
    Emm, Emv  = aux.get_energy_cg (args.m+"/geom_and_esurf.txt")
    print ("Emm = " + str(Emm) + ", Ems = " + str(Emv) + "...")
    df_M = pd.read_csv (args.m+"/energydump_1.mc", sep='\|', names=["energy", "mm_tot", "ms_tot", "timestep"], engine='python', skiprows=0)
    Nmm_M        = df_M["mm_tot"].values [-2000:]
    Nms_M        = df_M["ms_tot"].values [-2000:]
    net_energy_M = df_M["energy"].values [-2000:]
    net_energy_M = net_energy_M - np.mean( net_energy_M )

    U_CG   = Nmm_T*Emm + Nms_T*Emv
    U_T    = net_energy_T
    dU     = np.mean ( (U_CG - np.mean(U_CG)) - (U_T - np.mean(U_T)) )

    ############ BAR routine 
    # numerator: sample target configuration

    print ("<du/dlambda>_m - <du/dlambda>_t = ",np.mean (Nmm_M - Nmm_T))
    print ("N_mm (target) = " + str(np.mean (Nmm_T)) + ", N_mm (model) = " + str(np.mean(Nmm_M)) + ".")
