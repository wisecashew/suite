#!/usr/bin/env python3

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import argparse 
import aux
import re
import os 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide a potential energy surface.') 
parser.add_argument('--rg-file', dest='frg', metavar='RG_DATA', action='store', type=str, help='Name of RG data file.')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)
parser.add_argument('--rh-file', dest='frh', metavar='RH_DATA', action='store', type=str, help='Name of RH file.') 

args = parser.parse_args()   

if __name__=="__main__":
    
    # store the potential energy surface in a variable  
    U = args.U
    dop_list = aux.dir2dop (os.listdir(U))  
    print (dop_list[0:7]) 
    # go to the relevant potential energy, extract Rg value that corresponds with relevant temperature. 
    U_str  = "U = " + U + ":"
    rg_str = "Rg\^2:"
    rh_str = "invRh:"
    e_str  = "Error:"
    T_str  = "T:"
    
    
    fig = plt.figure()
    ax  = plt.axes()

    for N in dop_list[0:7]:

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
                print ("rg_list, hit")
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
        
        

        print ("T_list: ", T_list, flush=True)

        plt.plot ( T_list, np.sqrt(rg_list)*np.asarray(rh_list), marker='o', markeredgecolor='k' ) 
    

    plt.legend(["N = "+str(j) for j in dop_list[0:7] ] )
    ax.set_xscale ('log')
    ax.set_xlabel ( "Temperature (reduced)")
    ax.set_ylabel ( "$\sqrt{ \\langle R_g ^2 \\rangle }\cdot \\langle \\frac{1}{R_h} \\rangle$" ) 
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    ax.yaxis.get_major_locator().set_params(integer=True)
    
    plt.savefig ("rg_over_rh_"+"U_"+U+".png", dpi=800)
    
    if args.sp:
        plt.show()


