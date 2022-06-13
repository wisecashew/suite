#!/usr/bin/env python3

import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import argparse 
import aux
import re
import os 



''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', metavar='N', action='store', type=int, help='Provide a degree of polymerization.')
parser.add_argument('--rg-file', dest='frg', metavar='RG_DATA', action='store', type=str, help='Name of RG data file.')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)
parser.add_argument('--rh-file', dest='frh', metavar='RH_DATA', action='store', type=str, help='Name of RH file.') 

args = parser.parse_args()   

if __name__=="__main__":
    
    # store the potential energy surface in a variable  
    U_list = aux.dir2U ( os.listdir(".")) 
    N = args.dop 
    
    fig = plt.figure(1)
    ax  = plt.axes()
    
    rg_str = "Rg:"
    rh_str = "Rh:"
    e_str  = "Error:"
    T_str  = "T:"
    
    ax.set_xscale ( 'log' )
    ax.set_xlabel ( "Temperature (reduced)" )
    ax.set_ylabel ( "$\\langle R_g \\rangle / \\langle R_h \\rangle$" ) 
    # go to the relevant potential energy, extract Rg value that corresponds with relevant temperature. 
    
    i = 0  
    for U in U_list:
    
        U_str  = "U = " + U + ":"
        match_flag=False
        f = open (args.frg+"_"+str(N), 'r') 
        for line in f:
            if re.match ( U_str, line ):
                # print ("N is : " + str(N))
                # print ("line is: \"" + line[:-1] + "\"")
                match_flag = True
                continue
            
            if re.match( rg_str, line) and match_flag: 
                rg_list = line.strip().split() 
                rg_list = list ( map (float, rg_list[1:]) ) 
                continue
            
            elif re.match ( e_str, line) and match_flag:
                err_list = line.strip().split() 
                err_list = list ( map (float, err_list[1:]) ) 
                continue
            
            elif re.match ( T_str, line ) and match_flag: 
                T_list = line.strip().split() 
                T_list = list ( map (float, T_list[1:]) ) 
                break
        
        if not match_flag:
            print ("This potential energy was not found in N = {}.".format(N)) 
            quit()
    
        match_flag = False 
        f = open (args.frh + "_" + str(N), 'r' )
        for line in f:
            if re.match ( U_str, line ):
                # print ( "N is: " + str(N) ) 
                # print (" line is \" " + line[:-1] + "\"")
                match_flag = True 
                continue 
            if re.match ( rh_str, line) and match_flag: 
                rh_list = line.strip().split() 
                rh_list = list ( map (float, rh_list[1:]) ) 
                continue 

            elif re.match (e_str, line) and match_flag: 
                err_list = line.strip().split() 
                err_list = list ( map (float, err_list[1:] ) ) 
                continue 

            elif re.match ( T_str, line ) and match_flag: 
                T_list = line.strip().split() 
                T_list = list ( map (float, T_list[1:]) ) 
                break
        
        if not match_flag:
            print ("This potential energy was not foundin N = {}.".format(N) ) 
            quit()
        
        
        plt.plot ( T_list, np.asarray(rg_list)/np.asarray(rh_list), marker='o', markeredgecolor='k', color=cm.copper(i/9), label="_nolegend_", linestyle='-' ) 
        i += 1

    
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1))
    
    cbar = plt.colorbar (sm, orientation='vertical') 
    cbar.set_ticks ( [0,1] ) 
    cbar.set_ticklabels ( ["Weakest", "Strongest"] ) 
    cbar.ax.set_ylabel ( "Strength of aligned \nmonomer solvent interactions", fontsize=14, rotation=270 ) 
    plt.yticks (fontsize=12)
    plt.savefig ( "DOP_"+str(args.dop)+"_rg_over_rh.png", dpi=800)
    
    if args.sp:
        plt.show()


