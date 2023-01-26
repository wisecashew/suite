#!/usr/licensed/anaconda3/2020.7/bin/python

import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis")
import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    # U_list = aux.dir2U ( os.listdir(".") ) 
    U_list = ["U9"]
    PLOT_DICT = {}
    PLOT_DICT["U"]      = []
    PLOT_DICT["T"]      = []
    PLOT_DICT["U"]      = []
    PLOT_DICT["<E>"]    = []
    PLOT_DICT["Eps_mm"] = []
    PLOT_DICT["Eps_ms"] = []
    PLOT_DICT["Eps_ms/Eps_mm"] = []
    PLOT_DICT["MIX"]    = []
    PLOT_DICT["MIX/kT"] = []
    i=0
    ms_max = 25*2+(args.dop-2)*24
    for U in U_list:
        Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n = aux.get_energy (str(U)+"/geom_and_esurf.txt")
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True)
        energy_mean = np.asarray([])
        energy_list = np.asarray([])
        eps_ms_list = np.asarray([])
        eps_ms_err  = np.asarray([])
        eps_ms_mean = np.asarray([])
        eps_mm_list = np.asarray([])
        eps_mm_err  = np.asarray([])
        eps_mm_mean = np.asarray([])
        mix_mean    = np.asarray([])
        mix_list    = np.asarray([])
        beta_mix_mean = np.asarray([])
        beta_mix_list = np.asarray([])
        # temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        temperatures = [0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
        
        for temp in temperatures: 
            skip = 0
            energy_list = np.asarray   ([])
            eps_ms_list = np.asarray   ([])
            eps_mm_list = np.asarray   ([])
            mix_list    = np.asarray   ([])
            beta_mix_list = np.asarray ([])
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
                f = (df["mm_aligned" ].values[-2000:]*Emm_a  + df["mm_naligned" ].values[-2000:]*Emm_n )/(df["mm_tot"].values[-2000:])
                g = (df["ms1_aligned"].values[-2000:]*Ems1_a + df["ms1_naligned"].values[-2000:]*Ems1_n)/(df["ms1_tot"].values[-2000:])
                h = g - 0.5*f
                k = h/temp
                p = df["energy"].values[-2000:]
                eps_ms_list = np.hstack ( ( eps_ms_list, np.mean ( g ) ) )
                eps_mm_list = np.hstack ( ( eps_mm_list, np.mean ( f ) ) )
                mix_list    = np.hstack ( ( mix_list, np.mean ( h ) ) )
                beta_mix_list = np.hstack ( ( beta_mix_list, np.mean ( k ) ) )
                energy_list   = np.hstack ( ( energy_list, np.mean ( p ) ) )
            eps_ms_err  = np.hstack ( (eps_ms_err,  (np.std(eps_ms_list) / np.sqrt( len(num_list) ) ) ) )
            eps_ms_mean = np.hstack ( (eps_ms_mean, np.mean(eps_ms_list) ) )
            eps_mm_err  = np.hstack ( (eps_mm_err,  (np.std(eps_mm_list) / np.sqrt( len(num_list) ) ) ) )
            eps_mm_mean = np.hstack ( (eps_mm_mean, np.mean(eps_mm_list) ) )
            mix_mean    = np.hstack ( (mix_mean, np.mean ( mix_list ) ) )
            beta_mix_mean = np.hstack ( (beta_mix_mean, np.mean ( beta_mix_list ) ) )
            energy_mean   = np.hstack ( (energy_mean  , np.mean ( energy_list   ) ) )

        PLOT_DICT["U"].extend ( [U]*len (mix_mean) )
        PLOT_DICT["T"].extend ( temperatures )
        PLOT_DICT["Eps_mm"].extend (list(eps_mm_mean))
        PLOT_DICT["Eps_ms"].extend (list(eps_ms_mean))
        PLOT_DICT["Eps_ms/Eps_mm"].extend ( list(2*eps_ms_mean/eps_mm_mean) )
        PLOT_DICT["MIX"].extend (list(mix_mean))
        PLOT_DICT["MIX/kT"].extend (list(beta_mix_mean))
        PLOT_DICT["<E>"].extend (list(energy_mean))
        i += 1
        # mm_max.append( np.max(ms_list) ) 
        print("done!", flush=True)
    
    
    i=0
    chi_list = [0.1, 0.05, 0.01, 0.001, 0, -0.001, -0.01, -0.1, -0.2]
    rgba_color = cm.PiYG (divnorm(chi_list[i]))
    fig = plt.figure( num=1, figsize=(8,6) )
    ax  = plt.axes()
    ax.plot ( PLOT_DICT["T"], PLOT_DICT["MIX/kT"], linestyle='-', marker='o', markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
    ax.set_xscale ('log')
    ymin = 1
    df = pd.DataFrame.from_dict(PLOT_DICT)
    df.to_csv ("ENERGY-INFO.mc", index=False, sep='|')
    plt.savefig ("BETACHI-WHOLE.png", dpi=1200)
    fig = plt.figure (num=2, figsize = (8, 6))
    ax  = plt.axes()
    ax.plot ( PLOT_DICT["T"][8:], PLOT_DICT["MIX/kT"][8:], linestyle='-', marker='o', markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
    ax.set_xscale ('log')
    plt.savefig ("BETACHI-SUBSET.png", dpi=1200)
    fig = plt.figure (num=3, figsize = (8, 6))
    ax  = plt.axes()
    ax.plot ( PLOT_DICT["T"][8:], PLOT_DICT["Eps_ms/Eps_mm"][8:], linestyle='-', marker='o', markeredgecolor='k', linewidth=3, color=rgba_color, label="_nolabel_", markersize=10)
    ax.set_xscale ('log')
    plt.savefig ("RATIO.png", dpi=1200)
    fig = plt.figure(num=4, figsize=(8,6))
    ax = plt.axes()
    ax.plot ( PLOT_DICT["T"], PLOT_DICT["Eps_mm"], linestyle='-', marker='s', markeredgecolor='k', linewidth=3, color='steelblue', label="_nolabel_", markersize=10)
    ax.plot ( PLOT_DICT["T"], PLOT_DICT["Eps_ms"], linestyle='-', marker='^', markeredgecolor='k', linewidth=3, color='salmon', label="_nolabel_", markersize=10)
    ax.set_xscale ('log')
    plt.savefig ( "ENERGY-COMP.png", dpi=1200)
