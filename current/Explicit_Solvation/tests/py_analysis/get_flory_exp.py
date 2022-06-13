#!/usr/bin/env python3

import numpy as np 
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt 
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
parser.add_argument('-T', dest='T', action='store', type=float, help='Provide a temperature to find flory exponent.') 
# parser.add_argument('--excl-vol', dest='ev', action='store_true', help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--dump-file', dest='e', metavar='RG_DATA_X', action='store', type=str, help='Name of RG data file.')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)

args = parser.parse_args()   

if __name__=="__main__":
    
    # store the potential energy surface in a variable  
    U = args.U
    T = args.T
    
    print ("U = ", U)
    print ("T = ", T) 
    
    N_list = aux.dir2RGDATA ( os.listdir(".") ) 
    # print ( rgdatfiles )
    
    # go to the relevant potential energy, extract Rg value that corresponds with relevant temperature. 
    U_str  = "U = " + U + ":"
    rg_str = "Rg:"
    e_str  = "Error:"
    T_str  = "T:"
    rg_master_list = [] 
    for N in N_list[0:9]:
        match_flag=False
        f = open (args.e+"_"+str(N), 'r') 
        for line in f:
            if re.match ( U_str, line ):
                print ("N is : " + str(N))
                print ("line is: \"" + line[:-1] + "\"")
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
        
        idx = T_list.index(T) 
        rg_master_list.append (rg_list[idx]) 

    # outside of N loop 
    z = np.polyfit ( list(np.log (N_list[:9])), list(np.log(rg_master_list[0:9])), 1 ) 
    plt.figure(figsize=(8,6))
    plt.xscale('log')
    
    plt.plot (N_list[0:9], rg_master_list[0:9], marker='o', markeredgecolor='k')
    plt.plot (N_list[0:9], np.exp( z[0]*np.log(N_list[0:9]) + z[1]), marker='^', markeredgecolor='k', alpha=0.5)
    
    print ("Flory exponent is {}.".format(z[0]))

    plt.legend(["Simulation", "Fit"])
    plt.savefig("flory_plot_"+str(U)+"_"+str(T)+".png", dpi=800) 
    if args.sp:
        plt.show()


