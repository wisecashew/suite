#!/home/satyend/.conda/envs/data_analysis/bin/python

import pandas as pd 
import numpy as np 
import re
import copy
import os
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import multiprocessing 
import itertools 
import sys 
import copy
import time
import scipy.spatial.distance as ssd

sys.stdout.flush() 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

###########################################################################
###########################################################################
# Description: Looks at a string of form: "a | b | c | d"
# gets rid of the whitespace, and looks at all the other elements. 
# If it is an int, it holds it in loc 

def extract_loc_from_string(a_string):
    loc = [int(word) for word in a_string.split() if word.isdigit()]
    
    return loc     

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def get_frac (topology):
    
    f = open (topology)
    frac = "frac"
    

    for line in f:
        if re.findall (frac, line):
            r = re.findall("[0-9]+\.[0-9]+", line)
            break 

    f.close()
    
    return float(r[0])

def get_chi_fh (topology):
    
    f = open (topology)
    E_mm = "Emm_a"
    E_ms = "Ems_a"

    for line in f:
        if re.findall (E_mm, line):
            demixing = float( line[8:] )
        elif re.findall (E_ms, line):
            mixing   = float( line[8:] ) 

    chi = mixing - 0.5*demixing
    f.close()
    return chi 

def get_chi_entropy (topology):

    f = open (topology)
    Emm_a = "Emm_a"
    Emm_n = "Emm_n"
    Ems_a = "Ems_a"
    Ems_n = "Ems_n"

    for line in f:
        if re.findall (Emm_a, line):
            mm_a = float( line[8:] ) 
        elif re.findall ( Emm_n, line):
            mm_n = float( line[8:] ) 
        elif re.findall ( Ems_a, line):
            ms_a = float( line[8:] ) 
        elif re.findall ( Ems_n, line):
            ms_n = float( line[8:] ) 

    chi_a = ms_a - 0.5*mm_a
    chi_n = ms_n - 0.5*mm_n 
    f.close()

    return (chi_a, chi_n) 

def get_chi_cosolvent (topology):
    f = open (topology) 
    Emm_a  = "Emm_a"
    Emm_n  = "Emm_n"
    Ems1_a = "Ems1_a" 
    Ems2_a = "Ems2_a"
    frac   = "frac" 
    for line in f:
        if re.findall ( Emm_a, line):
            r = re.findall( "-[0-9]+\.[0-9]+|[0-9]+\.[0-9]+|-[0-9]+|[0-9]+", line)
            mm_a = float( r[0] )
        elif re.findall (Ems1_a, line):
            r = re.findall( "-[0-9]+\.[0-9]+|[0-9]+\.[0-9]+|-[0-9]+|[0-9]+", line)
            ms1_a = float( r[1] ) 
        elif re.findall ( Ems2_a, line):
            r = re.findall( "-[0-9]+\.[0-9]+|[0-9]+\.[0-9]+|-[0-9]+|[0-9]+", line)
            ms2_a = float ( r[1] ) 
        elif re.findall ( Emm_n, line):
            r = re.findall( "-[0-9]+\.[0-9]+|[0-9]+\.[0-9]+|-[0-9]+|[0-9]+", line)
            mm_n = float ( r[0] ) 
        elif re.findall ( frac, line): 
            r = re.findall( "-[0-9]+\.[0-9]+|[0-9]+\.[0-9]+|-[0-9]+|[0-9]+", line)
            frac_float = float ( r[0] ) 

    f.close()
    chi_1 = ms1_a - 0.5*mm_a 
    chi_2 = ms2_a - 0.5*mm_a 
    chi_3 = ms1_a - 0.5*mm_n

    return (chi_1, chi_2, chi_3, frac_float)

# End of function
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
###########################################################################
###########################################################################
# Description: Looks at a divident and divisor. It is a play on the original modulo function. 
# The regular % n return numbers from 0 to n-1.
# modified_modulo returns numbers from -n/2+1 to n/2

def modified_modulo(divident, divisor):
    midway = divisor/2
    if (divident%divisor > midway):
        result = (divident%divisor)-divisor 
        return result
    else:
        return divident%divisor         

# End of function 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: Custom function for my particular set of simulations. 
# it will give edge length of a simulation cell based of N i.e. the degree of 
# polymerization. 

def edge_length(N):
    
    adj = np.ceil(5*(N**(0.6)))
    if adj>N+2:
        return N+2
    else:
        return adj

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: When inside a directory, this function will figure out how many simulations
# were run. 

def dir2nsim (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        r = re.findall ("^energydump_[0-9]+.mc$", dir_name) 
        if ( r ):
            x = re.findall("\d+", r[0] ) 
            l.append ( int (x[0] ) )
    l.sort() 
    return l

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

############################################################################
############################################################################
# Description: When inside a directory, this function will find all the RG_DATA files
# and sort them according to degree of polymerization 

def dir2RGDATA (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        r = re.findall ("^RG_DATA_", dir_name)
        if (r):
            l.append ( int(dir_name[8:] ) )
    l.sort()
    return l

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def dir2set (list_of_dirs):
    l = []
    for dir_name in list_of_dirs:
        r = re.findall ("^set\d+",dir_name)
        if (r):
            l.append ( int(dir_name[3:] ) )
    l.sort()
    return l


############################################################################
############################################################################
# Description: When inside a directory, this function will figure out how
# many temperatures were sampled. 

def dir2dop (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        r = re.findall ("DOP_[0-9]+$", dir_name) 
        if ( r ):
            l.append ( int (r[0][4:] ) )
    l.sort() 
    return l

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: The polymer in the coordinate file is subjected to periodic boundary conditions
# so the coordinates cannot be processed directly. They need to be "unfucked" (unwrapped).

def unfuck_polymer(polymer, x, y, z): 
    unfucked_polymer = copy.copy(polymer) # np.asarray([polymer[0,:]])
     
    for i in range ( polymer.shape[0]-1 ) :
        diff = polymer[i+1,:] - polymer[i,:]
        
        for j in range(3):
            diff[j] = diff[j] if np.abs(diff[j])==1 else -1*np.sign(diff[j])
        
        unfucked_polymer[i+1] = unfucked_polymer[i]+diff # np.vstack( (unfucked_polymer, unfucked_polymer[i]+diff ) )
    
    return unfucked_polymer         


def unfuck_polymer_modified (polymer, x, y, z): 
    unfucked_polymer = copy.copy(polymer) # np.asarray([polymer[0,:]])
    oris = []
    # print ("polymer = ",polymer)
    for i in range ( polymer.shape[0]-1 ) :
        diff = polymer[i+1,:] - polymer[i,:]
        oris.append (int(polymer[i][3]))
        for j in range(3):
            diff[j] = diff[j] if np.abs(diff[j])==1 else -1*np.sign(diff[j])
        
        unfucked_polymer[i+1] = unfucked_polymer[i]+diff # np.vstack( (unfucked_polymer, unfucked_polymer[i]+diff ) )
    oris.append (int(polymer[polymer.shape[0]-1][3]))
    # print (oris)
    # quit()
    return (unfucked_polymer, oris)
# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: If the directory names are all floats, this function will get a list of floats from it.
# primarily used to get the range of temperatures. 

def dir2float (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        try:
            l.append(float(dir_name)) 
        except ValueError:
            continue
    l.sort()
    return l 

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: Given a directory with files labelled UX, I want to get a list of all of these U files
# and then have them sorted. 

def dir2U (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        if (re.match("^U\d+$", dir_name)):
            l.append(dir_name)
    
    l.sort()  
    l = sorted(l, key=lambda x: int(x[1:]) )
    return l

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: Given the dictionary which contains the coordinates of a polymer at all time points 
# in a trajectory, this function gives me the MEAN radius of gyration FOR THAT TRAJECTORY. 

def get_Rg (master_dict, xlen, ylen, zlen):
    
    N = master_dict [ next(iter(master_dict)) ][0].shape[0] 
    
    count = 0
    rg    = 0

    for key in master_dict: 
        coord_arr = unfuck_polymer ( master_dict[key][0], xlen, ylen, zlen )
        r_com = np.mean( coord_arr, axis=0) # get center of mass 
        offset = coord_arr - r_com 
        rg += np.sqrt(np.sum ( np.square (offset) )/ N) # added the np.sqrt
        count += 1

    return rg/count

#############################################################################
#############################################################################

def get_Rg2 (master_dict, xlen, ylen, zlen):
    
    N = master_dict [ next(iter(master_dict)) ][0].shape[0] 
    
    count = 0
    rg    = 0

    for key in master_dict: 
        coord_arr = unfuck_polymer ( master_dict[key][0], xlen, ylen, zlen )
        r_com = np.mean( coord_arr, axis=0) # get center of mass 
        offset = coord_arr - r_com 
        rg += np.sum (np.square (offset))/N # added the np.sqrt
        count += 1

    return rg/count


# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def get_Rg_list(master_dict, xlen, ylen, zlen):
    
    N = master_dict [ next(iter(master_dict)) ][0].shape[0] 
    
    rg_list    = []

    for key in master_dict: 
        coord_arr = unfuck_polymer ( master_dict[key][0], xlen, ylen, zlen )
        r_com = np.mean( coord_arr, axis=0) # get center of mass 
        offset = coord_arr - r_com 
        rg_list.append(  np.sum ( np.square (offset) )/ N  )
        
    return rg_list


#############################################################################
#############################################################################
# Description: get_rg_components 

def get_Rg_components ( coords_arr, xlen, ylen, zlen): 

    coords_arr = unfuck_polymer ( coords_arr, xlen, ylen, zlen ) 
    r_com = np.mean ( coords_arr, axis=0 ) 
    N = coords_arr.shape[0] 
    
    rsumx = np.sum( (coords_arr[:,0] - r_com[0])**2 )/N 
    rsumy = np.sum( (coords_arr[:,1] - r_com[1])**2 )/N
    rsumz = np.sum( (coords_arr[:,2] - r_com[2])**2 )/N
    
    comp_delta = ((rsumx**2) * (rsumy**2) + (rsumy**2) * (rsumz**2) + (rsumx**2) * (rsumz**2))/(rsumx**2+rsumy**2+rsumz**2)**2

    return comp_delta

def get_shape_param ( master_dict, xlen, ylen, zlen ):
    
    rdelta = [] 
    for key in master_dict:
        rdelta.append ( get_Rg_components ( master_dict[key][0], xlen, ylen, zlen ) ) 

    avg_comp = np.mean ( rdelta ) 
    return 1-3*avg_comp 


#############################################################################
#############################################################################
# Description: Given the dictionary which contains the coordinates of a polymer at all time points
# in a trajectory, this function gives me the hydrodynamic radius of that polymer FOR THAT TRAJECTORY. 

def get_Rh(master_dict, xlen, ylen, zlen):
    
    N = master_dict[ next(iter(master_dict)) ][0].shape[0]
    
    count = 0
    # Rh    = 0 
    inv_Rh = 0
    for key in master_dict: 
        coord_arr = unfuck_polymer(master_dict[key][0], xlen, ylen, zlen)
        count += 1 
        # Rh += 1/ ( np.sum( 1 / ssd.pdist ( coord_arr, 'euclidean' ) )/ (N*(N-1)/2) ) 
        inv_Rh += np.sum ( 1 / ssd.pdist (coord_arr, 'euclidean') ) / (N*(N-1)/2)  

    avg_inv_hydrodynamic_radius = inv_Rh/count  

    return avg_inv_hydrodynamic_radius 

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# Description: Given a trajectory, grab the coordinates of the polymer and load them into a dictionary. 
# this will allow easy access to these coordinates. 

def get_pdict(filename, starting_step, dop, x, y, z):
    
    st_b_str = "Dumping coordinates at step" 
    pmer_num_str = "Dumping coordinates of Polymer"
    start_str = "START"
    end_str_1 = "END" 
    end_str_2 = "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n" 
    
    starting_bool = False 
    # master_dict will be the dictionary which will contain polymer coordinates at each step     
    master_dict = {} 
    xlen, ylen, zlen = x, y, z
    
    step_flag     = 0
    pmer_flag     = 0
    end_step_flag = 0
    step_num      = 0

    # given a string, it will extract all numbers out in chronological order 
    # and put them in a numpy array 
    m_index = 0 
    f = open (filename, 'r')
    for line in f:
        if ( re.search(st_b_str, line) ):
            
            step_num = int ( ( extract_loc_from_string ( line.replace('.', ' ') ) )[0] )
            
            if ( step_num == starting_step ) or starting_bool :
                starting_bool = True
            else: 
                continue 
            # print ("in first if step_num is " + str(step_num) )
            
            master_dict [step_num] = {}
            
            step_flag     = 1
            pmer_flag     = -1
            end_step_flag = 0
            continue 

        elif ( re.search(start_str, line) ) and starting_bool:
            # print ("inside elif1, step_num is " + str(step_num) )
            continue

        elif ( re.search(pmer_num_str, line) ) and starting_bool:
            pmer_flag += 1
            master_dict[step_num][pmer_flag] = np.empty ( (dop,3) ) 
            continue 

        elif ( re.search(end_str_1, line) ) and starting_bool: 
            end_step_flag = 1
            step_flag     = 0
            pmer_flag     = -1 
            continue 

        elif ( re.search(end_str_2, line) ) and starting_bool:
            m_index = 0 
            # print ("index completed: ", step_num)
            continue 

        elif ( starting_bool ):
            # print ( "step_num is " + str(step_num) )
            monomer_coords                   = extract_loc_from_string ( line ) 
            master_dict[step_num][pmer_flag][m_index] = monomer_coords[0:-1] 
            m_index += 1
            continue
    
    f.close() 
    return master_dict

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def get_pdict_modified (filename, starting_step, dop, x, y, z):
    
    st_b_str = "Dumping coordinates at step" 
    pmer_num_str = "Dumping coordinates of Polymer"
    start_str = "START"
    end_str_1 = "END" 
    end_str_2 = "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n" 
    
    starting_bool = False 
    # master_dict will be the dictionary which will contain polymer coordinates at each step     
    master_dict = {} 
    xlen, ylen, zlen = x, y, z
    
    step_flag     = 0
    pmer_flag     = 0
    end_step_flag = 0
    step_num      = 0

    # given a string, it will extract all numbers out in chronological order 
    # and put them in a numpy array 
    m_index = 0 
    f = open (filename, 'r')
    for line in f:
        if ( re.search(st_b_str, line) ):
            
            step_num = int ( ( extract_loc_from_string ( line.replace('.', ' ') ) )[0] )
            
            if ( step_num == starting_step ) or starting_bool :
                starting_bool = True
            else: 
                continue 
            # print ("in first if step_num is " + str(step_num) )
            
            master_dict [step_num] = {}
            
            step_flag     = 1
            pmer_flag     = -1
            end_step_flag = 0
            continue 

        elif ( re.search(start_str, line) ) and starting_bool:
            # print ("inside elif1, step_num is " + str(step_num) )
            continue

        elif ( re.search(pmer_num_str, line) ) and starting_bool:
            pmer_flag += 1
            master_dict[step_num][pmer_flag] = np.empty ( (dop,4) ) 
            continue 

        elif ( re.search(end_str_1, line) ) and starting_bool: 
            end_step_flag = 1
            step_flag     = 0
            pmer_flag     = -1 
            continue 

        elif ( re.search(end_str_2, line) ) and starting_bool:
            m_index = 0 
            # print ("index completed: ", step_num)
            continue 

        elif ( starting_bool ):
            # print ( "step_num is " + str(step_num) )
            monomer_coords                   = extract_loc_from_string ( line ) 
            master_dict[step_num][pmer_flag][m_index] = monomer_coords 
            m_index += 1
            continue
    
    f.close() 
    return master_dict

#############################################################################
#############################################################################
# Description: This function is meant for the parallelism of plot_rg_parallel. 
# Given a potential energy surface, temperature, degree of polymerization, name of coordinates file 
# and the starting index, this function will go in and get the mean RADIUS OF GYRATION for THAT trajectory. 

def infiltrate_coords_get_rg ( U, T, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num)+".mc" 
    edge = edge_length (dop)
    # if T < 2.0:
    #     starting_index = 55000000
    # elif T == 2.5 or T == 5.0:
    #     starting_index = 90000000

    master_dict = get_pdict (filename, starting_index, dop, edge, edge, edge)
    # for key in master_dict:
    #    print ("key is ", key)
    #    print ("coordinates are: \n", master_dict[key][0] )
    # print ("Inside " + filename)
    rg = get_Rg(master_dict, edge, edge, edge) 
    return rg 

# End of function.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#############################################################################
#############################################################################
# Description: This function is meant for the parallelism of plot_rh_parallel. 
# Given a potential energy surface, temperature, degree of polymerization, name of coordinates file
# and the starting index, this function will go in and get the HYDRODYNAMIC RADIUS for THAT trajectory. 

def infiltrate_coords_get_rh ( U, T, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num) + ".mc" 
    # print (filename, flush=True)
    edge = edge_length (dop)
    master_dict = get_pdict (filename, starting_index, dop, edge, edge, edge) 
    rh = get_Rh(master_dict, edge, edge, edge) 

    return rh 

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#############################################################################
#############################################################################
# ALL GET FLORY FUNCTIONS NEED TO BE TESTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DATE: MAY 26, 2022

def get_flory(U, T, dop_list):
    
    rg_mean = [] 

    for dop in dop_list:
        filename = U+"/DOP_"+str(dop)+"/"+T+"/coords.txt" 
        master_dict = get_pdict (filename, dop+2, dop+2, dop+2)

        rg    = np.asarray([]) 
        steps = np.asarray([]) 
        
        for key in master_dict: 
            rg = np.hstack ( (rg, get_Rg2 (master_dict[key][0], dop+2, dop+2, dop+2) ) ) 
            steps = np.hstack ( steps, key )
        rg_mean = np.hstack ( rg_mean, np.mean(rg[args.s:] ) )

    rg_mean = np.array (rg_mean)
    
    L = np.polyfit ( np.log(dop_list), np.log(rg_mean), 1, w = np.sqrt(dop_list) ) 
    
    return L[0]/2

# NOT YET TESTED@!!!!!!!!!!!!!!!!!!!!!!!!!!!
#############################################################################
#############################################################################
# ALL GET FLORY FUNCTIONS NEED TO BE TESTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DATE: MAY 26, 2022

def plot_flory_multi ( U, T, starting_index, coords_files, show_plot_bool ):

    fig = plt.figure ( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16) 
    i=0
    Tmax = [] 

    rg_mean = np.array([])
    rg_std  = np.array([]) 

    dop_list = dir2nsim( os.listdir ( U ) )
    num_list = [1, 2, 3, 4, 5]
    coords_files = "coords"

    pool = multiprocessing.Pool ( processes=len(dop_list) )
    with pool as p:
        results = p.starmap ( infiltrate_coords_get_rg, zip (itertools.repeat(U), itertools.repeat(T), \
            dop_list, num_list, itertools.repeat(coords_files, itertools.repeat(starting_index) ) )  )

    pool.close()
    pool.join() 

    rg = np.asarray([])
    for r in results:
        rg = np.hstack (rg, np.mean (r) )
    rg_mean = np.hstack ( ( rg_mean, np.mean (rg) ) )
    # rg_std  = np.hstack ( ( rg_std, np.std(rg)/np.sqrt( len(num_list) ) ) )

    L = np.polyfit ( np.log(dop_list), np.log(rg_mean), 1, w = np.sqrt(dop_list) ) 

    return L[0]/2 

# NOT YET TESTED@!!!!!!!!!!!!!!!!!!!!!!!










##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION over
# a. a single U values 
# b. a single T value 
# c. In a single plot 
# d. For a given degree of polymerization 
# MEANT FOR TESTING, primarily. Need to check if things are working as expected 
# at a smaller, more testable scale. 
# This function has been parallelized to improve wall-clock time!

def plot_rg_parallelized_singletrajectory ( U, T, dop, starting_index, traj_num, coords_file ): 
    
    fig = plt.figure( figsize=(8,6) )
    ax = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16) 

    i    = 0
    
    edge = edge_length(dop)

    print ("Inside U = " + str(U) + ", and N = " + str(dop), flush=True) 
    filename = U+"/DOP_"+str(dop)+"/"+str(T)+"/"+coords_file+"_"+str(traj_num)+".mc"  
    master_dict = get_pdict ( filename, starting_index, dop, edge, edge, edge )

    rg_mean = get_Rg_list ( master_dict, edge, edge, edge )

    # for testing purposes: 
    
    ax.plot ( range(0, len(rg_mean)), rg_mean, marker='o', markeredgecolor='k', \
            markeredgewidth=1.5, linestyle='-', linewidth=0.5, label='_no_lengend') 
    ax.set_xscale('log')
    print ("mean rg = ", np.mean(rg_mean))
    # plt.savefig( "DOP_"+str(dop)+U+"_rg.png", dpi=800 ) 
    # plt.show() 

    return None

# End of function. 
#/////////////////////////////////////////////////////////////////////////


##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION over
# a. a single U values 
# b. range of T values 
# c. In a single plot 
# d. For a given degree of polymerization 
# MEANT FOR TESTING, primarily. Need to check if things are working as expected 
# at a smaller, more testable scale. 
# This function has been parallelized to improve wall-clock time!

def plot_rg_parallelized_single_U_dop (U, dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):
    
    U_list = [U]
        
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 


    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("RG_DATA_"+str(dop)+"_SINGULAR", "w") 

    for U in U_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop) + "...", flush=True )
        rg_mean = [] 
        rg_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            # pool = multiprocessing.Pool ( processes=len(mtemp_list[uidx]) ) # len(num_list)) 
            #with pool as p:
            results = pool_list[ shitty_dict[uidx] ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            # print (len(results))
            # pool.join() 

            print ("Pool has done its job. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rg_dict[mtemp_list[uidx][k]].append( results[k] )
        
            for T in np.unique (mtemp_list[uidx]):
                rg_mean.append( np.mean ( rg_dict[T] ) ) 
                rg_std.append ( np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) 
         
        print (rg_std)
        ax.errorbar   ( temperatures, rg_mean/np.sqrt(dop), yerr=rg_std/np.sqrt(dop), fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, color='darkgreen', label='_nolegend_' ) 
        
        f.write("Rg: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    
    f.close() 
    
    # plot Uexcl...
    if excl_vol_bool:
        temperatures = dir2float ( os.listdir( "Uexcl" +"/DOP_"+str(dop) ) )
        edge = edge_length (dop) 
        rg_mean = []
        rg_std  = [] 
        for T in temperatures:
            rg_list = [] 
            
            filename = "Uexcl/DOP_" + str(dop) + "/" + str(T) + "/" +coords_files + "_1"
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge ) 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
                r_com     = np.mean ( coord_arr, axis=0 ) 
                offset    = coord_arr - r_com
                rg_list.append ( np.sqrt ( np.sum ( np.square (offset) ) / dop ) ) 
            
            rg_mean.append ( np.mean (rg_list)/np.sqrt(dop) ) 
            rg_std.append  ( np.std  (rg_list)/np.sqrt(dop) ) 
        
        ax.errorbar ( temperatures, np.asarray(rg_mean), yerr=np.asarray(rg_std)/np.sqrt(5), fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1 )


    ax.set_xscale('log')
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g \\rangle/\sqrt{N}$", fontsize=18)     
    ax.set_xticks ( np.arange(np.min(Tmax), np.max(Tmax)+1, 1 ) )
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.savefig   ( "DOP_"+str(dop)+"_SINGULAR_rg.png", dpi=800)
    
    if show_plot_bool:
        plt.show() 

    return None

# End of function. 
#/////////////////////////////////////////////////////////////////////////


##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION over
# a. a single U values 
# b. range of T values 
# c. In a single plot 
# d. For a range of degrees of polymerization 
# MEANT FOR TESTING, primarily. Need to check if things are working as expected 
# at a smaller, more testable scale. 
# This function has been parallelized to improve wall-clock time!

def plot_rg_parallelized_single_U_multiple_dop (U, starting_index, excl_vol_bool, coords_files, show_plot_bool ):
    
    # U_list = [U]
        
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 
    
    dop_list = dir2dop ( os.listdir (U) ) 
    legend_list = [ "N = " + str(d) for d in dop_list ]

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("RG_DATA_"+str(dop)+"_SINGULAR", "w") 

    for dop in dop_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop) + "...", flush=True )
        rg_mean = [] 
        rg_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            num_list = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            results = pool_list[ shitty_dict[uidx] ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            print ("Pool has done its job. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rg_dict[mtemp_list[uidx][k]].append( results[k] )
        
            for T in np.unique (mtemp_list[uidx]):
                rg_mean.append( np.mean ( rg_dict[T] ) ) 
                rg_std.append ( np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) 
         
        ax.errorbar   ( temperatures, rg_mean/np.sqrt(dop), yerr=rg_std/np.sqrt(dop), fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, color='darkgreen', label='_nolegend_' ) 
        
        f.write("Rg: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    
    f.close() 
    
    # plot Uexcl...
    if excl_vol_bool:
        temperatures = dir2float ( os.listdir( "Uexcl" + "/DOP_" + str(dop) ) )
        edge = edge_length (dop) 
        rg_mean = []
        rg_std  = [] 
        
        for T in temperatures:
            rg_list = []     
            filename = "Uexcl/DOP_" + str(dop) + "/" + str(T) + "/" +coords_files + "_1"
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge ) 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
                r_com     = np.mean ( coord_arr, axis=0 ) 
                offset    = coord_arr - r_com
                rg_list.append ( np.sqrt ( np.sum ( np.square (offset) ) / dop ) ) 
            
            rg_mean.append ( np.mean (rg_list)/np.sqrt(dop) ) 
            rg_std.append  ( np.std  (rg_list)/np.sqrt(dop) ) 
        
        ax.errorbar ( temperatures, np.asarray(rg_mean), yerr=np.asarray(rg_std)/np.sqrt(5), fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1 )


    ax.set_xscale('log')
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g \\rangle/\sqrt{N}$", fontsize=18)     
    ax.set_xticks ( np.arange(np.min(Tmax), np.max(Tmax)+1, 1 ) )
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.savefig   ( "DOP_"+str(dop)+"_SINGULAR_rg.png", dpi=800)
    
    if show_plot_bool:
        plt.show() 

    return None

# End of function. 
#/////////////////////////////////////////////////////////////////////////



##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION over 
# a. Multiple U values 
# b. Multiple T values 
# c. In a single plot
# d. For a given degree of polymerization values
# Important and relatively involved/complex. 
# This function has been parallelized to improve wall-clock time! 

def plot_fh_rg_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):
    
    U_list = dir2U ( os.listdir (".") )
    PLOT_DICT = {} 
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    rg_max = 0 
    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("RG_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rg_mean = [] 
        rg_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {} 
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            # pool = multiprocessing.Pool ( processes=len(mtemp_list[uidx]) ) # len(num_list)) 
            #with pool as p:
            results = pool_list[ shitty_dict[uidx] ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            # print (len(results))
            # pool.join() 

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rg_dict[mtemp_list[uidx][k]].append( results[k] )
        
            for T in np.unique (mtemp_list[uidx]):
                rg_mean.append( np.mean ( rg_dict[T] ) ) 
                rg_std.append ( np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) 
        
        
        # ax.errorbar   ( temperatures, np.asarray(rg_mean)/dop, yerr=np.asarray(rg_std)/dop, fmt='o', markeredgecolor='k', \
        #             linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        
        PLOT_DICT[U] = (np.asarray(rg_mean), np.asarray(rg_std) ) 

        if rg_max < np.max (rg_mean):
            rg_max = np.max (rg_mean) 

        f.write("Rg^2: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 

    i=0
    for U in U_list:
        ax.errorbar ( temperatures, PLOT_DICT[U][0]/rg_max, yerr=PLOT_DICT[U][1]/rg_max, fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1

    # plot Uexcl...
    if excl_vol_bool:
        f = open ("RG_DATA_"+str(dop), 'a')
        f.write ( "U = Uexcl:\n" ) 
        print ("Inside U = Uexcl and N = " + str(dop), flush=True )
        temperatures_excl = dir2float ( os.listdir( "Uexcl" +"/DOP_"+str(dop) ) )
        edge = edge_length (dop) 
        rg_mean = []
        rg_std  = []
        for T in temperatures_excl:
            rg_list = [] 
            filename = "Uexcl/DOP_" + str(dop) + "/" + str(T) + "/" +coords_files 
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge ) 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
                r_com     = np.mean ( coord_arr, axis=0 ) 
                offset    = coord_arr - r_com
                rg_list.append ( np.sum ( np.square (offset)/dop ) ) 
           
            # print( rg_list ) 
            rg_mean.append( np.mean (rg_list) )
            rg_std.append( np.std  (rg_list) ) 
        # print ("rg_mean is ", rg_mean )
        # print ("rg_max is ", rg_max )
        ax.errorbar ( temperatures, np.ones(len(temperatures))*rg_mean[0]/rg_max, yerr=0, fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1 )
        ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
        
        f.write ("Rg^2: ")
        for j in range (len(temperatures) ):
            f.write ("{:.2f} ".format(rg_mean[0])) 
        f.write("\n")
        f.write("Error: ")
        for j in range(len (temperatures) ):
            f.write ( "{:.2f} ". format(0))
        f.write("\n")
        f.write("T: ") 
        for elem in temperatures:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.flush()
        f.close() 
    ########################################

    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Poorest", "Best"] )
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ("Quality of solvent", fontsize=18, rotation=270)
    ax.set_xscale('log')
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g^2 \\rangle/ \\langle R_g ^2 \\rangle _{\\mathrm{max}}$", fontsize=18)     
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    ax.set_yticks ( np.linspace(0,1,11) )
    plt.savefig   ( "DOP_"+str(dop)+"_multiple_rg.png", dpi=1000)
    
    if show_plot_bool:
        plt.show() 

    return None

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION over 
# a. Multiple U values 
# b. Multiple T values 
# c. In a single plot
# d. For a given degree of polymerization values
# Important and relatively involved/complex. 
# This function has been parallelized to improve wall-clock time! 

def plot_entropy_rg_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):
    
    U_list = dir2U ( os.listdir (".") )
    PLOT_DICT = {} 
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    rg_max = 0 
    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("RG_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rg_mean = [] 
        rg_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            # pool = multiprocessing.Pool ( processes=len(mtemp_list[uidx]) ) # len(num_list)) 
            #with pool as p:
            results = pool_list[ shitty_dict[uidx] ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            # print (len(results))
            # pool.join() 

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rg_dict[mtemp_list[uidx][k]].append( results[k] )
        
            for T in np.unique (mtemp_list[uidx]):
                rg_mean.append( np.mean ( rg_dict[T] ) ) 
                rg_std.append ( np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) 
        
        if rg_max < np.max (rg_mean):
            rg_max = np.max(rg_mean) 

        PLOT_DICT [U] = (np.asarray(rg_mean), np.asarray(rg_std))
        
        f.write("Rg^2: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 
    
    i=0
    for U in U_list:
        
        ax.errorbar   ( temperatures, PLOT_DICT[U][0]/rg_max, yerr=PLOT_DICT[U][1]/rg_max, fmt='o', markeredgecolor='k', \
                    linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i += 1

    # plot Uexcl...
    if excl_vol_bool:
        f = open ("RG_DATA_"+str(dop), 'a')
        f.write ( "U = Uexcl:\n")
        temperatures_excl = dir2float ( os.listdir( "Uexcl" +"/DOP_"+str(dop) ) )
        edge = edge_length (dop) 
        rg_mean = []
        rg_std  = [] 
        for T in temperatures_excl:
            rg_list = [] 
            
            filename = "Uexcl/DOP_" + str(dop) + "/" + str(T) + "/" +coords_files
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge ) 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
                r_com     = np.mean ( coord_arr, axis=0 ) 
                offset    = coord_arr - r_com
                rg_list.append ( np.sum ( np.square (offset)/dop ) ) 
            
            rg_mean.append ( np.mean (rg_list) ) 
            rg_std.append  ( np.std  (rg_list) ) 
        
        ax.errorbar ( temperatures, np.ones(len(temperatures))*rg_mean[0]/(rg_max), yerr=0 , fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1 )
        ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
        f.write ("Rg^2: ")
        for j in range (len (temperatures) ):
            f.write ( "{:.2f} ".format (rg_mean[0] ) ) 
        f.write("\n") 
        f.write("Error: ") 
        for j in range(len(temperatures) ):
            f.write( "{:.2f} ".format(0) ) 
        f.write ("\n") 
        for elem in temperatures:
            f.write ( "{:.2f} ".format(elem) ) 
        f.write ("\n") 
        f.flush()
        f.close() 
    ########################################

    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Weakest", "Strongest"] ) 
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ("Strength of aligned \nmonomer-solvent interactions", fontsize=18, rotation=270)
    ax.set_xscale('log')
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g^2 \\rangle/ \\langle R_g ^2 \\rangle _{\\mathrm{max}}$", fontsize=18)     
    ax.set_yticks (np.linspace(0, 1, 11)) 
    plt.savefig   ( "DOP_"+str(dop)+"_multiple_rg.png", dpi=1000)
    
    if show_plot_bool:
        plt.show() 

    return None

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#############################################################################
############################################################################


def plot_entropy_scaled_rg_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):
    
    U_list = dir2U ( os.listdir (".") )

    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    
    # plot Uexcl...
    if excl_vol_bool:
        temperatures = dir2float ( os.listdir( "U1" +"/DOP_"+str(dop) ) )
        temperatures_excl = dir2float ( os.listdir( "Uexcl" +"/DOP_"+str(dop) ) )
        edge = edge_length (dop) 
        rg_mean = []
        rg_std  = [] 
        for T in temperatures_excl:
            rg_list = [] 
            
            filename = "Uexcl/DOP_" + str(dop) + "/" + str(T) + "/" +coords_files
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge ) 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
                r_com     = np.mean ( coord_arr, axis=0 ) 
                offset    = coord_arr - r_com
                rg_list.append ( np.sum ( np.square (offset) ) / dop ) 
            
            rg_mean.append ( np.mean (rg_list)/dop ) 
            rg_std.append  ( np.std  (rg_list)/dop ) 
        
        ax.errorbar ( temperatures, np.ones(len(temperatures)), yerr=rg_std[0]/(rg_mean[0]*np.sqrt(5))*np.ones(len(temperatures)) , fmt='^', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1 )
        ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
        
        rg0 = rg_mean[0]
    
    ##########################################################


    i = 0 
    Tmax = [] 


    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("SCALED_RG_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rg_mean = [] 
        rg_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            # pool = multiprocessing.Pool ( processes=len(mtemp_list[uidx]) ) # len(num_list)) 
            #with pool as p:
            results = pool_list[ shitty_dict[uidx] ] .starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            # print (len(results))
            # pool.join() 

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rg_dict[mtemp_list[uidx][k]].append( results[k] )
        
            for T in np.unique (mtemp_list[uidx]):
                rg_mean.append( np.mean ( rg_dict[T] ) ) 
                rg_std.append ( np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) 
        
        
        ax.errorbar   ( temperatures, np.asarray(rg_mean)/(rg0*dop), yerr=np.asarray(rg_std)/(rg0*dop), fmt='o', markeredgecolor='k', \
                    linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        
        f.write("Rg^2: ") 
        for elem in rg_mean: 
            f.write ( "{:.2f} ".format(elem))
        f.write ("\n") 
        f.write ("Error: ")
        for elem in rg_std: 
            f.write ( "{:.2f} ".format(elem) )
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures: 
            f.write ( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()  

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 

    
    ########################################

    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Weakest", "Strongest"] ) 
    cbar.ax.set_ylabel ("Strength of aligned \nmonomer-solvent interactions", fontsize=18, rotation=270)
    ax.set_xscale('log')
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g^2 \\rangle /  \\langle R_g^2 \\rangle _0 $", fontsize=18)     
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    plt.savefig   ( "DOP_"+str(dop)+"_scaled_multiple_rg.png", dpi=800)
    
    if show_plot_bool:
        plt.show() 

    return None

# End of function. 


##########################################################################
##########################################################################
# Description: This function is meant for plotting HYDRODYNAMIC RADIUS over 
# a. Multiple U values 
# b. Multiple T values 
# c. In a single plot
# d. For a given degree of polymerization values
# Important and relatively involved/complex. 
# This function has been parallelized to improve wall-clock time! 

def plot_rh_parallelized ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
        
    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("INV_RH_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write("U = " + U + ":\n")
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rh_mean = [] 
        rh_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rh_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list( np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) )  ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rh_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1}

        for uidx in range(3):
            
            results = pool_list[ shitty_dict[uidx] ] .starmap ( \
                    infiltrate_coords_get_rh, zip( itertools.repeat(U), mtemp_list[uidx],\
                    mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rh_dict[mtemp_list[uidx][k]].append(results[k]) 
        
            for T in np.unique (mtemp_list[uidx]):
                rh_mean.append( np.mean ( rh_dict[T] ) ) 
                rh_std.append ( np.std  ( rh_dict[T] ) / np.sqrt( ntraj_dict[T] ) ) 

        # ax.errorbar   ( temperatures, rh_mean*np.sqrt(dop), yerr=rh_std, fmt='o', markeredgecolor='k', \
        #        linestyle='-', elinewidth=1, capsize=0, linewidth=1, \
        #        color=cm.copper(i/9), label='_nolegend_' ) 
        
        f.write("invRh: ") 
        for elem in rh_mean:
            f.write ("{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("Error: ") 
        for elem in rh_std:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.flush() 
        i+=1 
    

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 
    
    # plot Uexcl...
    if excl_vol_bool:
        f = open ("INV_RH_DATA_"+str(dop), 'a')
        f.write ("U = Uexcl:\n")
        temperatures_excl = dir2float ( os.listdir( "Uexcl/DOP_"+str(dop) ) )
        edge = edge_length (dop)
        rh_list = []
        for T in temperatures_excl:
            filename = "Uexcl/DOP_"+str(dop)+"/"+str(T)+"/"+coords_files
            edge = edge_length (dop) 
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge )
            
            count = 0 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge)
                rh_list.append ( np.sum ( 1/ ssd.pdist( coord_arr, 'euclidean' ) )/(dop*(dop-1)/2)  )
            
            rh_mean = np.mean (rh_list) 
            rh_std  = np.std  (rh_list) 

        ax.errorbar ( temperatures, np.ones (len(temperatures))*rh_mean*np.sqrt(dop), yerr=np.sqrt(dop)*np.ones (len(temperatures))*rh_std/np.sqrt(20), fmt='^', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1 ) 
        plt.legend  ( ["Athermal solvent"], loc='best', fontsize=12 ) 
        
        f.write("invRh: ")
        for j in range(len(temperatures)):
            f.write ("{:.2f} ".format(rh_mean) ) 
        f.write("\n")
        f.write("Error :")
        for j in range(len(temperatures)):
            f.write ("{:.2f} ".format(0))
        f.write("\n")
        for elem in temperatures:
            f.write( "{:.2f} ".format(elem)) 
        f.write("\n")
        f.flush()
        f.close()
####################################
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )

    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Poor solvent", "Good solvent"] ) 
    cbar.ax.set_ylabel ("Solvent quality", fontsize=18, rotation=270)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle \\frac{1}{R_h} \\rangle \cdot \sqrt{N}$", fontsize=18)     
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    plt.savefig   ( "DOP_"+str(dop)+"_inv_rh.png", dpi=800)
    if show_plot_bool:
        plt.show() 
    
    return None

#############################################################################
#############################################################################




##########################################################################
##########################################################################
# Description: This function is meant for plotting HYDRODYNAMIC RADIUS over 
# a. Multiple U values 
# b. Multiple T values 
# c. In a single plot
# d. For a given degree of polymerization values
# Important and relatively involved/complex. 
# This function has been parallelized to improve wall-clock time! 

def plot_entropy_rh_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
        
    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("INV_RH_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write("U = " + U + ":\n")
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rh_mean = [] 
        rh_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rh_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list( np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) )  ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rh_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1}

        for uidx in range(3):
            
            results = pool_list[ shitty_dict[uidx] ] .starmap ( \
                    infiltrate_coords_get_rh, zip( itertools.repeat(U), mtemp_list[uidx],\
                    mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rh_dict[mtemp_list[uidx][k]].append(results[k]) 
        
            for T in np.unique (mtemp_list[uidx]):
                rh_mean.append( np.mean ( rh_dict[T] ) ) 
                rh_std.append ( np.std  ( rh_dict[T] ) / np.sqrt( ntraj_dict[T] ) ) 

        ax.errorbar   ( temperatures, np.asarray(rh_mean)*np.sqrt(dop), yerr=np.asarray(rh_std)*np.sqrt(dop), fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, \
                color=cm.copper(i/9), label='_nolegend_' ) 
        
        f.write("invRh: ") 
        for elem in rh_mean:
            f.write ("{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("Error: ") 
        for elem in rh_std:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.flush() 
        i+=1 
    

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 
    
    # plot Uexcl...
    if excl_vol_bool:
        f = open ("INV_RH_DATA_"+str(dop), 'a') 
        f.write ("U = Uexcl:\n")
        temperatures_excl = dir2float ( os.listdir( "Uexcl/DOP_"+str(dop) ) )
        edge = edge_length (dop)
        rh_list = []
        for T in temperatures_excl:
            filename = "Uexcl/DOP_"+str(dop)+"/"+str(T)+"/"+coords_files
            edge = edge_length (dop) 
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge )
            
            count = 0 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge)
                rh_list.append (  np.sum ( 1/ ssd.pdist( coord_arr, 'euclidean' ) )/(dop*(dop-1)/2)  )
            
            rh_mean = np.mean (rh_list) 
            rh_std  = np.std  (rh_list) 

        ax.errorbar ( temperatures, np.ones (len(temperatures))*rh_mean*np.sqrt(dop), yerr=np.ones (len(temperatures))*rh_std/(np.sqrt(dop)*np.sqrt(20)), fmt='^', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1 ) 
        plt.legend  ( ["Athermal solvent"], loc='best', fontsize=12 ) 

        f.write("invRh: ")
        for j in range(len(temperatures)):
            f.write("{:.2f} ".format(rh_mean) )
        f.write("\n")
        f.write("Error: ")
        for j in range(len(temperatures)):
            f.write  ("{:.2f} ".format(0)) 
        f.write("\n")
        f.write("T: ")
        for elem in temperatures:
            f.write ("{:.2f} ".format(elem))
        f.write ("\n")
        f.flush()
        f.close()


######################################
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )

    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Weakest", "Strongest"] ) 
    cbar.ax.tick_params (labelsize=14)
    cbar.ax.set_ylabel ("Strength of aligned \nmonomer-solvent interactions", fontsize=18, rotation=270)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\left\\langle \\frac{1}{R_h} \\right\\rangle \cdot \sqrt{N}$", fontsize=18)     
    ax.set_xscale ("log")
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
    # ax.set_xticks ( temperatures )
    plt.savefig   ( "DOP_"+str(dop)+"_rh.png", dpi=800)
    if show_plot_bool:
        plt.show() 
    
    return None

#############################################################################
#############################################################################



def plot_fh_rh_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
        
    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )# len(num_list)) 
    pool2 = multiprocessing.Pool ( processes=5 )
    
    pool_list = [pool1, pool2]
    
    f = open("RH_DATA_"+str(dop), "w") 

    for U in U_list:
        f.write("U = " + U + ":\n")
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        rh_mean = [] 
        rh_std  = [] 
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rh_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list( np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) )  ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rh_dict[T] = []


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1}

        for uidx in range(3):
            
            results = pool_list[ shitty_dict[uidx] ] .starmap ( \
                    infiltrate_coords_get_rh, zip( itertools.repeat(U), mtemp_list[uidx],\
                    mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                rh_dict[mtemp_list[uidx][k]].append(results[k]) 
        
            for T in np.unique (mtemp_list[uidx]):
                rh_mean.append( np.mean ( rh_dict[T] ) ) 
                rh_std.append ( np.std  ( rh_dict[T] ) / np.sqrt( ntraj_dict[T] ) ) 

        ax.errorbar   ( temperatures, np.asarray(rh_mean)*np.sqrt(dop), yerr=np.asarray(rh_std)*np.sqrt(dop), fmt='o', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1, \
                color=cm.copper(i/9), label='_nolegend_' ) 
        
        f.write("invRh: ") 
        for elem in rh_mean:
            f.write ("{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("Error: ") 
        for elem in rh_std:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.write("T: ") 
        for elem in temperatures:
            f.write( "{:.2f} ".format(elem) ) 
        f.write("\n") 
        f.flush() 
        i+=1 
    

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join() 

    f.close() 
    
    # plot Uexcl...
    if excl_vol_bool:
        f = open("INV_RH_DATA_"+str(dop), 'a')
        f.write ("U = Uexcl:\n")
        temperatures_excl = dir2float ( os.listdir( "Uexcl/DOP_"+str(dop) ) )
        edge = edge_length (dop)
        rh_list = []
        for T in temperatures_excl:
            filename = "Uexcl/DOP_"+str(dop)+"/"+str(T)+"/"+coords_files
            edge = edge_length (dop) 
            master_dict = get_pdict ( filename, 0, dop, edge, edge, edge )
            
            count = 0 
            for key in master_dict:
                coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge)
                rh_list.append (  np.sum ( 1/ ssd.pdist( coord_arr, 'euclidean' ) )/(dop*(dop-1)/2)  )
            
            rh_mean = np.mean (rh_list) 
            rh_std  = np.std  (rh_list) 

        ax.errorbar ( temperatures, np.ones (len(temperatures))*rh_mean*np.sqrt(dop), yerr=np.ones (len(temperatures))*rh_std/(np.sqrt(dop)*np.sqrt(20)), fmt='^', markeredgecolor='k', \
                linestyle='-', elinewidth=1, capsize=0, linewidth=1 ) 
        plt.legend  ( ["Athermal solvent"], loc='best', fontsize=12 ) 
        
        f.write("invRh: ")
        for j in range(len(temperatures)):
            f.write ("{:.2f} ".format(rh_mean) )
        f.write("\n")
        for j in range(len(temperatures)):
            f.write ("{:.2f} ".format(0))
        f.write ("\n") 
        for elem in temperatures:
            f.write ("{:.2f} ".format(elem))
        f.write("\n")
        f.flush()
        f.close() 
######################################
    my_cmap = cm.copper 
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )

    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Poorest", "Best"] ) 
    cbar.ax.tick_params (labelsize=14)
    cbar.ax.set_ylabel ("Quality of solvent", fontsize=18, rotation=270)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\left\\langle \\frac{1}{R_h} \\right\\rangle \cdot \sqrt{N}$", fontsize=18)     
    ax.set_xscale ("log")
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
    # ax.set_xticks ( temperatures )
    plt.savefig   ( "DOP_"+str(dop)+"_rh.png", dpi=1000)
    if show_plot_bool:
        plt.show() 
    
    return None

#############################################################################
#############################################################################







##########################################################################
##########################################################################
# Description: This function is meant for plotting RADIUS OF GYRATION given 
# a. a certain U value 
# b. a certain T value 
# c. in a single plot
# d. for a given degree of polymerization values
# e. for a given trajectory
# # Information obtained from trajectory 

def print_rg_singletraj ( U, T, dop, traj_num, starting_index, coords_files ):

    rg = infiltrate_coords_get_rg ( U, T, traj_num, dop, coords_files, starting_index)
    
    print ("File name provided is: " + U + "/DOP_" + str(dop) + "/" + str(T) + "/" + coords_files + "_" + str(traj_num) ) 
    print ("Starting index is ", starting_index)
    print ("Mean radius of gyration is: ", rg)

    return None 
    
# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


##########################################################################
##########################################################################
# Description: This function is meant for plotting shape factor 
# a. given a certain U value 
# b. given a certain T value 
# c. given a certain N values 
# d. given a traj file 
# e. given a traj number
# print out the shape factor 

def shape_factor ( U, T, num, dop, coords_file, starting_index ):
    
    filename = U+"/DOP_"+str(dop)+"/"+str(T)+"/"+coords_file+"_"+str(num)
    edge     = edge_length (dop)
    master_dict  = get_pdict( filename, starting_index, dop, edge, edge, edge) 
    count = 0 
    shape_term = 0
    for key in master_dict:
        coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
        r_com = np.mean ( coord_arr, axis=0 )
        offset = coord_arr - r_com 
        rgx = np.sqrt ( np.sum( np.square (offset)[:,0]   )/dop )
        rgy = np.sqrt ( np.sum( np.square (offset)[:,1]   )/dop ) 
        rgz = np.sqrt ( np.sum( np.square (offset)[:,2]   )/dop ) 
        shape_term += ( (rgx**2)*(rgy**2) + (rgy**2)*(rgz**2) + (rgx**2)*(rgz**2) )/( (rgx**2) + (rgy**2) + (rgz**2) )**2 
        count += 1 

    delta = 1 - 3*shape_term/count 

    return delta 
    
##########################################################################
##########################################################################

def single_sim_flory_exp ( U, T, num, dop, coords_file, starting_index, delta ):
    
	filename = U+"/DOP_"+str(dop)+"/"+str(T)+"/"+coords_file+"_"+str(num)+".mc"
	edge     = edge_length (dop)
	master_dict  = get_pdict( filename, starting_index, dop, edge, edge, edge) 
	offset_list = []

	for key in master_dict:
		coord_arr    = unfuck_polymer ( master_dict[key][0], edge, edge, edge ) 
		delta_coords = coord_arr [0:dop-delta] - coord_arr [delta:]
		offset = list(np.linalg.norm ( delta_coords, axis=1 ) **2 )
		offset_list.extend(offset) 

	return np.mean (offset_list) 

##########################################################################
##########################################################################

def get_e2e(coord_arr, xlen, ylen, zlen): 
    coord_arr = unfuck_polymer(coord_arr, xlen, ylen, zlen) 
    
    e2e = coord_arr[-1] - coord_arr[0] 
    return e2e 

# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#############################################################################

def get_e2e_data ( U, T, dop, traj_num, starting_index, coords_files, show_plot ): 

    import tidynamics 
    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/" + coords_files + "_" + str(traj_num)
    xlen = edge_length(dop)
    ylen = edge_length(dop)
    zlen = edge_length(dop)
    
    start = time.time() 
    print ("Making the dict...", flush=True)  
    master_dict = get_pdict ( filename, starting_index, dop, xlen, ylen, zlen )  
    print ("Made the dict!", flush=True) 
    end = time.time() 
    print ("Time to make dict = {:.2f}".format(end-start) )
    e2e = [] 
    print ("Making e2e array...", flush=True)
    for key in master_dict: 
        e2e.append ( get_e2e (master_dict[key][0], xlen, ylen, zlen ) )
    
    print ("Made the array!", flush=True) 
    end2 = time.time() 
    print ("Time to make e2e array = {:.2f}".format(end2-end) )
    # auto_corr = [] 
    
    print ("Making the auto_correlation function...", flush=True)
    auto_corr = tidynamics.acf ( e2e ) 
    print ("Made the acf!", flush=True)
    end3 = time.time()
    print ("Time to make acf = {:.2f}".format(end3-end2) )
    auto_corr = auto_corr/auto_corr[0] 
    
    print ("auto_corr length is: ", len(auto_corr))

    delays = np.arange ( 0, int(len(auto_corr)/2))*1000
    plt.plot   ( delays, auto_corr[0:int(len(auto_corr)/2)], marker='o', markeredgecolor='k', markerfacecolor='darkgreen', linestyle='-' ) 
    plt.xlabel ( "$\delta$" ) 
    plt.ylabel ( "$ \\frac { \\langle R(\delta) R(0) \\rangle } { \\langle R(0) \\cdot R(0) \\rangle }$" ) 
    plt.savefig ( U+"T_"+str(T)+"DOP_"+str(dop) + ".png", dpi=1000) 
     
    if show_plot:
         plt.show() 
    
    return None 


# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#############################################################################
#############################################################################
# Description This function meant for plotting the shape parameter
# a. multiple U values 
# b. Multiple T values 
# c. In a single plot 
# d. For a given degree of polymerization 

def plot_fh_shape_parameter_parallelized_single_dop_all_U_all_T ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") ) 

    fig = plt.figure ( figsize=(8,6) ) 
    ax = plt.axes ( ) 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    i = 0
    Tmax = [] 

    # instantiating pool 
    pool1 = multiprocessing.Pool ( processes=50 ) 
    pool2 = multiprocessing.Pool ( processes=5 )

    pool_list = [pool1, pool2]

    for U in U_list:
        pass


    return None
# End of function. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#############################################################################
#############################################################################
# Description This function meant for obtaining the order parameter
# a. multiple U values 
# b. Multiple T values 
# c. In a single plot 
# d. For a given degree of polymerization 

def obtain_order_parameter ( U, N, T, ortn_file_name, idx, starting_index, norm ):
    
    Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
            3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
            6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
            9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
            12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
            15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
            18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }
    start_str = "START for Step"
    end_str   = "END"

    f = open(U+"/DOP_"+str(N)+"/"+str(T)+"/"+ortn_file_name+"_"+str(idx)+".mc", 'r')

    extract_orr = False
    start_bool  = False
    monomer_order  = []
    oparam_list    = []
    alignment      = [] 
    
    if T < 2:
        starting_index = 60000000

    for line in f:
        if re.match ( start_str, line ):
            a = re.search ("\d+", line)
            extract_orr = True
            if int ( a.group(0) ) == starting_index:
                start_bool = True
            count          = 0

        elif re.match ( end_str, line ) and start_bool:
            extract_orr = False
            oparam_list.append ( np.mean(monomer_order) )
            monomer_order.clear()
            alignment.clear() 

        elif extract_orr and start_bool:
            or_list    = extract_loc_from_string ( line ) [1:] # [1:]  these takes all the orientations of the solvent molecules
            monomer_or = extract_loc_from_string ( line ) [0]
            count = 0
            for cnum in or_list:
                alignment.append( np.dot (Or2Dir[monomer_or], Or2Dir[cnum]) )
                # count += 1
            if len(alignment) == 0:
                # monomer_order.append (0)
                continue
            else:
                # if norm == 'z':
                    # monomer_order.append (alignment/26)
                # elif norm == 'mm':
                    # monomer_order.append (alignment/count)
                monomer_order.append (  np.mean( np.asarray(alignment)**2  ) - (np.mean(alignment))**2 ) 

    f.close()
    fluctuation = np.mean(oparam_list) 
    return fluctuation 

##########################################################################
##########################################################################

def obtain_monomer_alignment ( U, N, T, coords_files, idx, starting_index, norm ):

    Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
            3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
            6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
            9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
            12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
            15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
            18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }

    filename = U + "/DOP_" + str (N) + "/" + str(T) + "/" + coords_files + "_" + str(idx) + ".mc"
    edge = edge_length (N)
    master_dict = get_pdict_modified ( filename, starting_index, N, edge, edge, edge )
    alignment_lists = []
    monomer_alignment = []
    monomer_monomer = 0
    count = 0
    
    for key in master_dict:
        (coord_arr, orientations) = unfuck_polymer_modified ( master_dict[key][0], edge, edge, edge)
        # print (master_dict[key])
        # print (len(orientations))
        # exit()
        monomer_alignment.clear() 
        monomer_monomer   = 0
        # identify all neighbors
        for i in range(N-1):
            for j in range(i+1, N):
                diff = coord_arr[j][:-1]-coord_arr[i][:-1]
                # print ("diff = ",diff)
                if (np.abs(int(diff[0])) == 1 or int(diff[0]) == 0) and (np.abs(int(diff[1])) == 1 or int(diff[1]) == 0) and (np.abs(int(diff[2])) == 1 or int(diff[2]) == 0):
                    monomer_alignment.append( np.dot (Or2Dir[orientations[i]], Or2Dir[orientations[j]]) )
                    # 
        if len(monomer_alignment) == 0: # monomer_monomer == 0:
            continue
        else:
            # if norm == 'z':
            #     alignment_lists.append ( monomer_alignment/26 )
            # elif norm == 'mm':
            #    alignment_lists.append ( monomer_alignment/monomer_monomer )
            alignment_lists.append (np.mean( np.asarray(monomer_alignment)**2) - (np.mean(monomer_alignment))**2  )

    fluctuations = np.mean (alignment_lists) 
    return fluctuations


###############################################################################
###############################################################################

def get_neighbors (polymer):
    v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
    [-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
    [1, 0, 1], [-1, 0, 1], [1, 0, -1], \
    [-1, 0, -1], [1, 1, 1], \
    [-1, 1, 1], [1, -1, 1], \
    [1, 1, -1], [-1, -1, 1], \
    [-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
    v = [np.asarray(u) for u in v] 
    v = np.asarray (v) 
    # print (polymer)
    fss_solvent_set = set()
    for p in polymer:
        for v_ in v:
            check = v_ + p
            if np.any(np.all(check == polymer, axis=1)):
                pass
            else: 
                fss_solvent_set.add(tuple( check ))
    sss_solvent_set = set() 
    for s in fss_solvent_set:
        for v_ in v:
            check = v_ + np.asarray (s) 
            if np.any(np.all(check==polymer, axis=1)) or (tuple(check) in fss_solvent_set):
                pass
            else:
                sss_solvent_set.add(tuple( check ))

    return list(fss_solvent_set), list(sss_solvent_set)
