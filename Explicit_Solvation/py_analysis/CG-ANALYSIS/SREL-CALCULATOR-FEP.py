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
    Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n  = aux.get_energy    ("TARGET/geom_and_esurf.txt")
    energies_t  = np.array ([Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n])
    df_T = pd.read_csv("TARGET/energydump_1.mc", sep='\|', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
    Nmm_T         = df_T["mm_tot" ].values[-2000:]
    Nms_T         = df_T["ms1_tot"].values[-2000:] + df_T["ms2_tot"].values[-2000:]
    net_energy_T  = df_T["energy" ].values[-2000:] 
    net_energy_T  = net_energy_T - np.mean( net_energy_T )
    energies_cg   = [-0.1, 0.0] # aux.get_energy_cg (args.m+"/geom_and_esurf.txt")

    #######################################
    for e_mm in np.arange(-0.3, 0.3, 0.05):# [0, -0.1, -0.15, -0.2, -0.25, -0.3]: 
        energies_cg[0] = e_mm
        print ("For e_mm = " + str(e_mm), end =', ')
        U_CG   = Nmm_T*energies_cg[0] + Nms_T*energies_cg[1]
        U_CG   = U_CG - np.mean(U_CG) 
        U_T    = net_energy_T
        dU     = np.mean ( U_CG - U_T )
        dfree_energy = -k*T*np.log ( np.mean ( np.exp ( -beta* (U_CG - U_T) ) ) )
        S_rel = dU - dfree_energy
        print ("S_rel = " + str(S_rel))
