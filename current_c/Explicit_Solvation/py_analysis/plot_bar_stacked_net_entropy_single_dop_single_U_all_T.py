#!/usr/licensed/anaconda3/2020.7/bin/python

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


parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide energy surface.') 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')

args = parser.parse_args()


if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = [args.U]
    fig, ax = plt.subplots(figsize=(8,6))
    # instantiate plt figure 
    # fig.figure( figsize=(8,6) )
    
    PLOT_DICT  = {}
    ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    # starting_index = args.s
    

    ms_max = 25*2+(args.dop-2)*24
    for U in U_list:
        
        temperatures = aux.dir2float (os.listdir(str(U)+"/DOP_32"))
        temperatures = [0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 1.0, 2.5, 10.0, 100.0]
        temperatures.sort()
        mma_list = np.asarray([])
        mma_mean = np.asarray([])
        mmn_list = np.asarray([])
        mmn_mean = np.asarray([])
        msa_list = np.asarray([])
        msa_mean = np.asarray([])
        msn_list = np.asarray([])
        msn_mean = np.asarray([])

        for T in temperatures: 
            
            mma_list = np.asarray ([]) 
            mmn_list = np.asarray ([])
            msa_list  = np.asarray ([]) 
            msn_list  = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
                mma_list  = np.hstack ( (mma_list, df["mm_aligned"  ].values[-2000:]*2/args.dop  ) )
                mmn_list  = np.hstack ( (mmn_list, df["mm_naligned" ].values[-2000:]*2/args.dop  ) )
                msa_list  = np.hstack ( (msa_list, df["ms1_aligned" ].values[-2000:]/args.dop ) )
                msn_list  = np.hstack ( (msn_list, df["ms1_naligned"].values[-2000:]/args.dop ) )

            msa_mean = np.hstack ( (msa_mean, np.mean(msa_list) ) )
            msn_mean = np.hstack ( (msn_mean, np.mean(msn_list) ) )
            mma_mean = np.hstack ( (mma_mean, np.mean(mma_list) ) )
            mmn_mean = np.hstack ( (mmn_mean, np.mean(mmn_list) ) )

        PLOT_DICT [U] = (mma_mean, mmn_mean, msa_mean, msn_mean)
        print("done!", flush=True)
    

    for U in U_list:
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][0], color ='darkred', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][1], bottom = PLOT_DICT[U][0], color='lightcoral', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][2], bottom = PLOT_DICT[U][0]+PLOT_DICT[U][1], color='steelblue', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][3], bottom = PLOT_DICT[U][0]+PLOT_DICT[U][1]+PLOT_DICT[U][2], color='lightskyblue', width=0.8, edgecolor='k')

    ###
    # write data points out for text later 
    g = open ("m1sn.dat", 'w')
    g.write ("m-m-a	m-m-n	m-s-a	m-s-n\n")
    for i in range(len (PLOT_DICT[U][0])):
        g.write ("{}	{}	{}	{}\n".format(PLOT_DICT[U][0][i], PLOT_DICT[U][1][i], PLOT_DICT[U][2][i], PLOT_DICT[U][3][i] ) )
    g.close()
    ax.legend(["$m$-$m$ (aligned)", "$m$-$m$ (misaligned)", "$m$-$s$ (aligned)", "$m$-$s$ (misaligned)"], ncol=2, loc="upper center", fancybox=True, bbox_to_anchor=(0.5, 1.13))
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=14, direction="in", left="off", labelleft="on", pad=5 )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    ax.axhline (y=0, c='k', linewidth=1)
    ax.minorticks_on()
    ax.set_ylim((-1 , 30))
    ax.set_xticks (np.arange(len(temperatures)))
    print ([str(i) for i in temperatures])
    ax.set_xticklabels ([str(i) for i in temperatures])
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0f}'))
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    
    plt.savefig   (args.pn + ".png", dpi=1200)

