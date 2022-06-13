#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np 
import re 
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib.cm as cm 
import os
import argparse 
import aux
import multiprocessing 


parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the average shape parameter for this trajectory.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-U', metavar='U', dest='U', type=str, action='store', help='enter a potential energy surface.') 
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number.', default=100)
parser.add_argument('--coords-file', dest='e', action='store', type=str, metavar='coords', help='Name of coords file to parse through.') 
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on a screen.', default=False)

args = parser.parse_args() 

if __name__=="__main__":
    
    U_list = [args.U] 

    # instantiate plt figure 
    plt.figure ( figsize=(8,6) ) 
    dop = args.dop 
    xlen = aux.edge_length (dop)
    ylen = aux.edge_length (dop)
    zlen = aux.edge_length (dop) 

    # instantiate some pertinent variables 
    i = 0
    
    for U in U_list: 
        
        print ( "Currently plotting out stuff in U = " + str(U) + "...", end = ' ', flush=True ) 
        temperatures = aux.dir2float ( os.listdir ( U+"/DOP_"+str(dop) ) ) 
        shape_parameter = [] 
        for T in temperatures: 
            num_list = aux.dir2nsim ( os.listdir ( U+"/DOP_"+str(dop)+"/"+str(T) ) ) 
            avg_shape = 0
            for num in num_list:
                filename = U+"/DOP_"+str(dop)+"/"+str(T)+"/"+args.e+"_"+str(num) 
                master_dict = aux.get_pdict (filename, args.s, args.dop, xlen, ylen, zlen) 
                avg_shape += aux.get_shape_param ( master_dict, xlen, ylen, zlen )
            shape_parameter.append ( avg_shape/len(num_list) ) 

        plt.plot ( temperatures, shape_parameter, marker='o', markeredgecolor='k', linestyle='-', color=cm.copper(i/9), label='_nolegend_' ) 
        i += 1
        print ("done!")

    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )

    ax = plt.axes() 

    plt.ylabel ( "Shape parameter ($\\delta ^*$)", fontsize=18)
    plt.xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_xscale ("log")
    plt.xticks ( temperatures) 
    cbar = plt.colorbar (sm, orientation='vertical') 
    cbar.set_ticks ( [0,1] ) 
    cbar.set_ticklabels ( ["Weakest", "Strongest"] ) 
    cbar.ax.set_ylabel ( "Strength of aligned \nmonomer-solvent interactions", fontsize=16, rotation=270 )
    plt.yticks (fontsize=12)
    plt.savefig ("DOP_"+str(dop)+"_shape_parameter.png", dpi=800)

    if args.sp:
        plt.show()

