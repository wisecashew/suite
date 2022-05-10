#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import re
import copy

###########################################################################
###########################################################################

def extract_loc_from_string(a_string):
    loc = [int(word) for word in a_string.split() if word.isdigit()]
    
    return np.asarray(loc)     

def modified_modulo(divident, divisor):
    midway = divisor/2
    if (divident%divisor > midway):
        result = (divident%divisor)-divisor 
        return result
    else:
        return divident%divisor         

#############################################################################
#############################################################################

def unfuck_polymer(polymer, x, y, z): 
    unfucked_polymer = np.asarray([polymer[0,:]])
    
    for i in range ( polymer.shape[0]-1 ) :
        diff = polymer[i+1,:] - polymer[i,:]
        
        for j in range(3):
            diff[j] = modified_modulo(diff[j], x)
        
        unfucked_polymer = np.vstack( (unfucked_polymer, unfucked_polymer[i]+diff ) )
    
    return unfucked_polymer         

#############################################################################
#############################################################################

def dir2float (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        try:
            l.append(float(dir_name)) 
        except ValueError:
            continue
    l.sort()
    return l 

#############################################################################
#############################################################################

def get_Rg(coord_arr, xlen, ylen, zlen):
    
    coord_arr = unfuck_polymer(coord_arr, xlen, ylen, zlen)
    
    r_com = np.mean(coord_arr, axis=0) 
    N = coord_arr.shape[0]
    rsum = 0
    
    

    for i in range(N): 
        rsum += np.linalg.norm( coord_arr[i,:]- r_com )**2 
    
    rsum = np.sqrt(rsum/N)
    
    return rsum

#############################################################################
#############################################################################

def get_Rh(master_dict, xlen, ylen, zlen, ignore_until):
    N = master_dict[0][0].shape[0] 
    
    inverse_distance = np.zeros(int(N*(N-1)/2)) 
    tot_steps = len(master_dict) 
    count = 0
    for key in master_dict:
        if (key/1000 < ignore_until):
            continue
        else: 
            coord_arr = unfuck_polymer(master_dict[key][0], xlen, ylen, zlen)
            k = 0
            count += 1 
            for i in range(N):
                for j in range(i+1,N):
                    inverse_distance[k] += 1/(np.linalg.norm( coord_arr[i] - coord_arr[j], 2 ) ) 
                    k += 1

    ens_average = np.sum(inverse_distance)/count  
    hydrodynamic_radius = 1/ ( ens_average/N**2 )  

    return hydrodynamic_radius 

#############################################################################
#############################################################################

def get_pdict(filename, x, y, z):
    f = open(filename, 'r')
    coord_file = f.readlines() 
    
    st_b_str = "Dumping coordinates at step" 
    pmer_num_str = "Dumping coordinates of Polymer"
    start_str = "START"
    end_str_1 = "END" 
    end_str_2 = "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n" 

    # master_dict will be the dictionary which will contain polymer coordinates at each step     
    master_dict = {} 
    xlen, ylen, zlen = x, y, z
    
    step_flag     = 0
    pmer_flag     = 0
    end_step_flag = 0
    step_num      = 0

    # given a string, it will extract all numbers out in chronological order 
    # and put them in a numpy array 
    
    for line in coord_file:
        if ( re.search(st_b_str, line)):
            step_num = int ( ( extract_loc_from_string ( line.replace('.', ' ') ) ) )
            master_dict [step_num] = {}
            
            step_flag     = 1
            pmer_flag     = -1
            end_step_flag = 0
            continue 

        elif ( re.search(start_str, line) ):
            continue

        elif ( re.search(pmer_num_str, line) ):
            pmer_flag += 1
            master_dict[step_num][pmer_flag] = np.empty ( (0,3) ) 
            continue 

        elif ( re.search(end_str_1, line) ): 
            end_step_flag = 1
            step_flag     = 0
            pmer_flag     = -1 
            continue 

        elif ( re.search(end_str_2, line) ):
            continue 

        else:
            monomer_coords                   = extract_loc_from_string ( line ) 
            master_dict[step_num][pmer_flag] = np.vstack ( (master_dict[step_num][pmer_flag], monomer_coords[0:-1] ) ) 
            continue
    
    return master_dict

#############################################################################
#############################################################################
