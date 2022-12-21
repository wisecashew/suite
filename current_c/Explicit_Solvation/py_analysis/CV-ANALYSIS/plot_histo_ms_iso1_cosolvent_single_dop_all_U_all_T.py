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
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')
parser.add_argument('-T', dest='T', metavar='temp', action='store', type=float, nargs='+', help='temperature ranges.')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm (0.001, vmin=-0.2, vmax=0.1 ) # this is for entropy 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    frac_list = []
    # U_list = ["U11"] 
    plt.figure( figsize=(8,6) )
    temperatures = args.T
    PLOT_DICT = {}
    i=0
    Tmax = []
    ms_max = 25*2+(args.dop-2)*24
    for temp in temperatures:
        print ("Currently plotting out stuff in T = " + str(temp) + "...", end=' ', flush=True)
        ms_list = np.asarray([])
        ms_err  = np.asarray([])
        ms_mean = np.asarray([])
        
        for U in U_list:
            frac_list.append (aux.get_frac(U+"/geom_and_esurf.txt"))
            skip = 0
            ms_list = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp) ) ) )
            for num in num_list: 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
                f = df["ms1_tot" ].values[-2000:]+df["ms2_tot"].values[-2000:] 
                ms_list = np.hstack ( ( ms_list, f ) )
                
            # ms_err  = np.hstack ( (ms_err,  (np.std(ms_list)/np.sqrt(30) ) ) )
            PLOT_DICT[U] = ms_list
            # print ("PLOT_DICT[" + str(temp) + "] = ",PLOT_DICT[temp])

        i += 1
        print("done!", flush=True)
    
    i=0
    ymin = 1
    ms_max = 1
    
    chi_list = [-0.01, -0.1, -0.2] # [0.1, 0.05, 0.01, 0.001, 0, -0.001, -0.01, -0.1, -0.2]
    for U in U_list:
        
        fig = plt.figure(i)
        ax = plt.axes()
        ax.minorticks_on()
        
        ax.hist ( PLOT_DICT[U], range=(400, 775), bins=np.arange(400, 775), density=True)
        ax.set_ylim (0, 0.1)
        fig.savefig (args.pn + "_" + str(frac_list[i]) + ".png", bbox_inches='tight', dpi=1200)
        # fig.clf ()
        i += 1


