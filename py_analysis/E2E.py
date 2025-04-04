#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np 
import re 
import matplotlib.pyplot as plt 
import pandas as pd
import aux

'''
This code will take in a trajectory file generated by my MonteCarlo engine and 
give you Radius of Gyration
'''

import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and end-to-end correlation function.")
parser.add_argument('-i', metavar=': input coordinate file from which Rg is calculated (coords.txt)', dest='i', action='store', help='enter address of coordinate file')
parser.add_argument('-e', metavar=': edge length of the box in which the simulation was conducted', type=int, dest='e', action='store', help='enter the edge length of the cubic simulation box')
args = parser.parse_args() 

def get_e2e(coord_arr, xlen, ylen, zlen): 
    coord_arr = aux.unfuck_polymer(coord_arr, xlen, ylen, zlen) 
    
    e2e = coord_arr[-1] - coord_arr[0] 
    return e2e 


if __name__ == "__main__":
    
    f = open( args.i , "r")
    # f = open("coords.txt", 'r')
    coord_file = f.readlines() 
    f.close() 
    
    st_b_str = "Dumping coordinates at step" 
    start_str = "START"
    pmer_num_str = "Dumping coordinates of Polymer"
    end_str_1 = "END" 
    end_str_2 = "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n" 
    
    # master_dict will be dictionary which will contain polymer coordinates at each step 
    # chunky data structure
    master_dict = {} 
    
    xlen, ylen, zlen = args.e, args.e, args.e
    
    step_flag = 0 
    pmer_flag = -1 
    end_step_flag = 0
    step_num = 0 
    
    # given a string, it will extract all numbers out in chronological order 
    # and put them in a numpy array 
    
    for line in coord_file: 
        if (re.search(st_b_str, line)):
            print (line)
            step_num = int(re.findall(r'\d+', line)[0])
            master_dict[step_num] = {}
            step_flag = 1
            pmer_flag = -1
            end_step_flag = 0
            continue
        
        elif (re.search(start_str, line)):
            continue

        elif (re.search(pmer_num_str, line)):
            pmer_flag += 1
            master_dict[step_num][pmer_flag] = np.empty ( (0,3) )
            
            continue
            
        elif (re.search(end_str_1, line)):
            end_step_flag = 1
            step_flag = 0 
            pmer_flag = -1            
            continue
        
        elif (re.search(end_str_2,line)):
            continue
        else:         
            monomer_coords = aux.extract_loc_from_string(line)
            # print(monomer_coords)
            master_dict[step_num][pmer_flag] = np.vstack( (master_dict[step_num][pmer_flag], monomer_coords[0:-1]) )
            continue
    
    e2e = [] 
    for key in master_dict: 
        e2e.append( get_e2e(master_dict[key][0], xlen, ylen, zlen ) ) 
    
    auto_corr = [] 
    delays = np.arange(0, int(len(e2e)/2))
    
    for i in range(int(len(e2e)/2)): 
        rsum=0 
        for j in range(len(e2e) - i ): 
            # print(j)
            rsum += np.dot(e2e[j], e2e[j+i]) 
        rsum = rsum/(len(e2e)-i) 
        # print(len(e2e)-i)
        auto_corr.append(rsum) 
    
    auto_corr = auto_corr/auto_corr[0] 
    
    # d = {'delays': delays, 'auto_corr': auto_corr}
    # pd.DataFrame(data=d)
    g = open("e2e_corr.dat", 'w')
    g.write("# delay    correlation\n") 
    for (d,c) in zip(delays, auto_corr):
        g.write("{}    {:.2f}\n".format(d, c))
    g.close()
        
    plt.plot(delays, auto_corr, '-')
    plt.xlabel("$\delta$")
    plt.ylabel("$ \\frac { \\langle R(\delta) R(0) \\rangle } { \\langle R(0) \\cdot R(0) \\rangle }$")
    plt.savefig("end2end_autocorr.png", dpi=1200) 
    # plt.show() 
     
