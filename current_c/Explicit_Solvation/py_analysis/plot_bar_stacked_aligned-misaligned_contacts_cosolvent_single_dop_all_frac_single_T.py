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
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

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
    starting_index = args.s
    

    ms_max = 25*2+(args.dop-2)*24
    for U in U_list:
        
        temperatures = aux.dir2float (os.listdir(str(U)+"/DOP_32"))
        temperatures = [0.01, 0.1, 0.5, 1.0,2.5, 5.0,  10.0, 50.0, 100.0]
        temperatures.sort()
        mma_list = np.asarray([])
        mma_err  = np.asarray([])
        mma_mean = np.asarray([])
        mmn_list = np.asarray([])
        mmn_err  = np.asarray([])
        mmn_mean = np.asarray([])
        msa_list = np.asarray([])
        msa_err  = np.asarray([])
        msa_mean = np.asarray([])
        msn_list = np.asarray([])
        msn_err  = np.asarray([])
        msn_mean = np.asarray([])

        for T in temperatures: 
            skip = args.s
            if T < 2.0:
                skip = 5500
            elif T == 2.5 or T == 5.0:
                skip = 9000
            mma_list = np.asarray ([]) 
            mmn_list = np.asarray ([])
            msa_list  = np.asarray ([]) 
            msn_list  = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=skip)
                mma_list  = np.hstack ( (mma_list, df["mm_aligned"].values*2/args.dop  ) )
                mmn_list  = np.hstack ( (mmn_list, df["mm_naligned"].values*2/args.dop ) )
                msa_list  = np.hstack ( (msa_list, df["ms1_aligned"].values/args.dop ) )
                msn_list  = np.hstack ( (msn_list, df["ms1_naligned" ].values/args.dop)) 
                
            msa_err  = np.hstack ( (msa_err , np.std (msa_list)/np.sqrt(30) ) )
            msa_mean = np.hstack ( (msa_mean, np.mean(msa_list) ) )
            msn_err  = np.hstack ( (msn_err , np.std (msn_list)/np.sqrt(30) ) )
            msn_mean = np.hstack ( (msn_mean, np.mean(msn_list) ) )
            mma_err  = np.hstack ( (mma_err , np.std (mma_list)/np.sqrt(30) ) )
            mma_mean = np.hstack ( (mma_mean, np.mean(mma_list) ) )
            mmn_err  = np.hstack ( (mmn_err , np.std (mmn_list)/np.sqrt(30) ) )
            mmn_mean = np.hstack ( (mmn_mean, np.mean(mmn_list) ) )
            

        PLOT_DICT [U] = (mma_mean, mmn_mean, msa_mean, msn_mean)
        ERROR_DICT[U] = (mma_err , mmn_err , msa_err, msn_err)
        # mm_max.append( np.max(ms_list) )
        print("done!", flush=True)
    

    for U in U_list:
        # plt.errorbar ( frac_list, PLOT_DICT[T]/ ms_max, yerr=ERROR_DICT[T]/ms_max, linewidth=1, fmt='none', ecolor=colors[i], capsize=2, color='k', label="_nolabel_")
        # plt.plot     ( frac_list, PLOT_DICT[T]/ ms_max, linestyle='-',  markeredgecolor='k', linewidth=3, color=colors[i], label="_nolabel_", markersize=10, marker='o')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][0], color='coral', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][1], bottom=PLOT_DICT[U][0], color='steelblue', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][2], bottom= PLOT_DICT[U][0]+PLOT_DICT[U][1], color='seagreen', width=0.8, edgecolor='k')
        ax.bar ( np.arange(len(temperatures)), PLOT_DICT[U][3], bottom= PLOT_DICT[U][0]+PLOT_DICT[U][1]+PLOT_DICT[U][2], color='oldlace', width=0.8, edgecolor='k')

    ###
    # write data points out for text later 
    g = open ("m1sn.dat", 'w')
    g.write ("m-ma	m-mn	m-sa	m-sn\n")
    for i in range(len (PLOT_DICT[U][0])):
        g.write ("{}	{}	{}	{}\n".format(PLOT_DICT[U][0][i], PLOT_DICT[U][1][i], PLOT_DICT[U][2][i], PLOT_DICT[U][3][i] ) )
    g.close()
    # my_cmap = cm.coolwarm
    # sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0) ) 
    # ax = plt.axes()
    ax.legend(["$m_1$-$m_1^a$", "$m_1$-$m_1 ^n$", "$m_1$-$s_1 ^a$", "$m_1$-$s_1 ^n$"], ncol=4)
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=14, direction="in", left="off", labelleft="on", pad=5 )
    ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
    # cbar = plt.colorbar(sm, orientation='vertical')
    # cbar.set_ticks ( [0, 0.25, 0.5, 0.75, 1.0] )
    # cbar.set_ticklabels( ["$10^{-2}$", "$10^{-1}$", "1", "$10^{1}$", "$10^2$"] ) 
    # cbar.set_ticks( [] )
    # cbar.ax.tick_params (labelsize=14)
    # cbar.set_ticklabels( [] )
    # cbar.ax.set_ylabel ( "Strength of better solvent \n", fontsize=16, rotation=270 ) 
    # ax.set_xscale('log')
    # ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    ax.axhline (y=0, c='k', linewidth=0.2)
    # ax.yaxis.get_major_locator().set_params(integer=True)
    ax.minorticks_on()
    # ax.set_xlim((-0.05, 1.05))
    ax.set_ylim((-1 , 30))
    ax.set_xticks (np.arange(len(temperatures)))
    print ([str(i) for i in temperatures])
    ax.set_xticklabels ([str(i) for i in temperatures])
    # ax.set_yticks (np.linspace (0.50,1,6))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0f}'))
    # plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    # ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    
    plt.savefig   (args.pn + ".png", bbox_inches='tight', dpi=1200)
