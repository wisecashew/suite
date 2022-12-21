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
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.')
parser.add_argument('-T', dest='T', action='store', type=float, nargs='+', help='Provide a temperature.')
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump') 
parser.add_argument('--png-name', dest='pn', metavar='png name', action='store', type=str, help='Name of image.', default='ms_plot')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False) 

args = parser.parse_args()

divnorm = matplotlib.colors.SymLogNorm ( 0.5, vmin=0, vmax=1.0 ) 

if __name__=="__main__":

    # get the entire list of potential energy surfaces 
    U_list = aux.dir2U ( os.listdir(".") ) 
    fig, ax = plt.subplots(figsize=(8,6))
    # instantiate plt figure 
    # fig.figure( figsize=(8,6) )
    
    PLOT_DICT  = {}
    ERROR_DICT = {}
    dop            = args.dop
    dump_file      = args.e
    starting_index = args.s
    
    print ("T = ", args.T)
    temperatures = [float(elem) for elem in args.T]
    temperatures.sort() 

    ms_max = 25*2+(args.dop-2)*24
    for T in temperatures:
        
        frac_list    = [] 
        print ("Currently plotting out stuff in T = " + str(T) + "...", end=' ', flush=True)
        ms1_list = np.asarray([])
        ms1_err  = np.asarray([])
        ms1_mean = np.asarray([])
        ms2_list = np.asarray([])
        ms2_err  = np.asarray([])
        ms2_mean = np.asarray([])
        mm_list  = np.asarray([])
        mm_err   = np.asarray([])
        mm_mean  = np.asarray([])

        for U in U_list: 
            frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
            ms1_list = np.asarray ([]) 
            ms2_list = np.asarray ([])
            mm_list  = np.asarray ([]) 
            num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U)+"/DOP_"+str(args.dop)+"/"+str(T) ) ) )

            for num in num_list: 
                # print (str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc") 
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(T)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=args.s) 
                ms1_list = np.hstack ( (ms1_list, df["ms1_tot"].values/args.dop ) )
                ms2_list = np.hstack ( (ms2_list, df["ms2_tot"].values/args.dop ) )
                mm_list  = np.hstack ( (mm_list , df["mm_tot" ].values*2/args.dop)) 
                
            ms1_err  = np.hstack ( (ms1_err , np.std (ms1_list)/np.sqrt(30) ) )
            ms1_mean = np.hstack ( (ms1_mean, np.mean(ms1_list) ) )
            ms2_err  = np.hstack ( (ms2_err , np.std (ms2_list)/np.sqrt(30) ) )
            ms2_mean = np.hstack ( (ms2_mean, np.mean(ms2_list) ) )
            mm_err   = np.hstack ( (mm_err  , np.std (mm_list)/np.sqrt(30)  ) ) 
            mm_mean  = np.hstack ( (mm_mean , np.mean(mm_list ) ) )

        PLOT_DICT [T] = (ms1_mean, ms2_mean, mm_mean)
        ERROR_DICT[T] = (ms1_err , ms2_err , mm_err )
        # mm_max.append( np.max(ms_list) )
        print("done!", flush=True)
    
    i=0
    if len(temperatures) == 1:
        colors = cm.coolwarm (np.linspace(0,1,3)) 
        i = 1
    else:
        colors = cm.coolwarm (np.linspace(0,1,len(temperatures)))

    for T in temperatures:
        # plt.errorbar ( frac_list, PLOT_DICT[T]/ ms_max, yerr=ERROR_DICT[T]/ms_max, linewidth=1, fmt='none', ecolor=colors[i], capsize=2, color='k', label="_nolabel_")
        # plt.plot     ( frac_list, PLOT_DICT[T]/ ms_max, linestyle='-',  markeredgecolor='k', linewidth=3, color=colors[i], label="_nolabel_", markersize=10, marker='o')
        ax.bar ( frac_list, PLOT_DICT[T][0], color='coral', width=0.05, edgecolor='k')
        ax.bar ( frac_list, PLOT_DICT[T][1], bottom=PLOT_DICT[T][0], color='steelblue', width=0.05, edgecolor='k')
        ax.bar ( frac_list, PLOT_DICT[T][2], bottom= PLOT_DICT[T][0]+PLOT_DICT[T][1], color='seagreen', width=0.05, edgecolor='k')
        i += 1

    ###
    # write data points out for text later 
    g = open ("m1sn.dat", 'w')
    g.write ("m-s1	m-s2	m-m	\n")
    for i in range(len (PLOT_DICT[T][0])):
        g.write ("{}	{}	{}	\n".format(PLOT_DICT[T][0][i], PLOT_DICT[T][1][i], PLOT_DICT[T][2][i] ) )
    g.close()
    # my_cmap = cm.coolwarm
    # sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0) ) 
    # ax = plt.axes()
    ax.legend(["$m_1$-$s_1$ contacts", "$m_1$-$s_2$ contacts", "$m_1$-$m_1$ contacts"], ncol=3)
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
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
    ax.set_xlim((-0.05, 1.05))
    ax.set_ylim((-1 , 30))
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # ax.set_yticks (np.linspace (0.50,1,6))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    
    plt.savefig   (args.pn + ".png", dpi=1000)
