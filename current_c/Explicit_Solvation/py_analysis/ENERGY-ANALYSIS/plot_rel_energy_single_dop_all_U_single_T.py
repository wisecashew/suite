#!/usr/licensed/anaconda3/2020.7/bin/python

import sys
sys.path.insert(0, "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis")
import pandas as pd 
import numpy as np 
import matplotlib
from matplotlib.ticker import StrMethodFormatter
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as tck
import argparse 
import aux 
import os 
import mpmath as mpm

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.')
parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    
    # instantiate plt figure 
    plt.figure( figsize=(8,6) )
    MASTER_DICT  = {} 
    MASTER_DICT["T"  ] = []
    MASTER_DICT["xc" ] = []
    MASTER_DICT["E"  ] = []
    MASTER_DICT["err_E"] = []
    MASTER_DICT["F"]     = []
    ENERGY_DICT  = {}
    ENERGY_ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    starting_index = args.s
    
    print ("T = ", args.T)
    temperatures = [float(elem) for elem in args.T]
    temperatures.sort() 

    Emm_a, Emm_n, Ems1_a, Ems1_n, Ems2_a, Ems2_n, Es1s2_a, Es1s2_n = aux.get_energy ("U1/geom_and_esurf.txt")
    ms_max = 25*2+(args.dop-2)*24
    for T in temperatures:
        
        frac_list    = [] 
        print ("Currently plotting out stuff in T = " + str(T) + "...", end=' ', flush=True)
        energy_list      = np.asarray ([])
        energy_err       = np.asarray ([])
        energy_mean      = np.asarray ([])
        free_energy_list = np.asarray ([])
        free_energy_err  = np.asarray ([])
        free_energy_mean = np.asarray ([])

        for U in U_list:
            frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
            energy_list = np.asarray ([]) 
            free_energy_list = np.asarray([])
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=args.s) 
                energy_list = np.hstack ( ( energy_list, df["energy"].values[-2000:] ) )
            energy_err  = np.hstack ( (energy_err , np.std (energy_list) / np.sqrt( len (num_list) ) ) )
            energy_mean = np.hstack ( (energy_mean, np.mean(energy_list) ) )
            free_energy_list = [mpm.exp(1/T * en_elem) for en_elem in energy_list] 
            free_energy_list = np.mean(free_energy_list)
            free_energy_mean = np.hstack ( (free_energy_mean, T*mpm.log(free_energy_list) ) )
            

        ENERGY_DICT [T]      = energy_mean
        ENERGY_ERROR_DICT[T] = energy_err
        MASTER_DICT ["T"].extend (list(np.ones(len(frac_list))*T))
        MASTER_DICT ["xc"].extend(frac_list) 
        MASTER_DICT ["E"].extend(list(energy_mean)) 
        MASTER_DICT ["err_E"].extend (list(energy_err))
        MASTER_DICT ["F"].extend (list(free_energy_mean))
        print("done!", flush=True)
    
    
    df = pd.DataFrame.from_dict (MASTER_DICT, orient='columns')
    # print (df)
    df.to_csv ("THERMODYNAMIC-INFO.csv", sep='|', index=False)
    
    
