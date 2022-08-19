#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import aux 
import re
import os 

parser = argparse.ArgumentParser(description="Get the plot of Rg versus temperatures from the RG_DATA file.")
parser.add_argument('-U', dest='U', action='store', type=str, help="Provide a potential energy surface.")
parser.add_argument('-T', dest='T', action='store', type=float, help='Provide a temperature to analyze.') 
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide a degree of polymerization.')
parser.add_argument('--dump-file', dest='e', metavar='RG_DATA_X', help="Provide the name of the dump file.")
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want an image to be rendered.') 

args = parser.parse_args() 

if __name__=="__main__":
    
    # instantiate variables
    U   = args.U
    T   = args.T 
    dop = args.dop

    # define regex'
    U_str  = "U = " + U + ":"
    rg_str = "Rg:"
    e_str  = "Error:"
    T_str  = "T:"

    rg_list = [] 

    match_flag = False
    f = open (args.e+"_"+str(args.dop), 'r')
    for line in f:
        if re.match (U_str, line):
            match_flag=True 
            continue 
        elif re.match ( rg_str, line ) and match_flag:
            rg_list = line.strip().split() 
            rg_list = list ( map (float, rg_list[1:]) ) 
            continue 
        elif re.match ( e_str, line ) and match_flag:
            err_list = line.strip().split()
            err_list = list ( map (float, err_list[1:]) ) 
            continue
        elif re.match ( T_str, line ) and match_flag:
            T_list = line.strip().split() 
            T_list = list ( map ( float, T_list[1:] ) ) 
            break
    
    print (T_list)
    plt.plot ( T_list, rg_list, marker='o', markeredgecolor='k', markerfacecolor='green', linestyle='-', linewidth=3 ) 
    plt.ylabel ( "$\\langle R_g \\rangle$" ) 
    plt.xlabel ( "Temperature (reduced)" ) 
    plt.xticks ( np.arange(0,11,1) ) 
    plt.legend ( ["U = " + str(U)] ) 
    plt.title ("N = " + str(args.dop) )
    plt.savefig ("rg_plot_dop"+str(dop)+"_"+str(U)+".png", dpi=800)
    if args.sp:
        plt.show()


