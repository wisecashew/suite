#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import re
import copy
import os
import matplotlib.pyplot as plt 
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import multiprocessing 
import itertools 


''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

###########################################################################
###########################################################################

def extract_loc_from_string(a_string):
    loc = [int(word) for word in a_string.split() if word.isdigit()]
    
    return np.asarray(loc)     

###########################################################################
###########################################################################

def modified_modulo(divident, divisor):
    midway = divisor/2
    if (divident%divisor > midway):
        result = (divident%divisor)-divisor 
        return result
    else:
        return divident%divisor         

#############################################################################
#############################################################################

def edge_length(N):
    
    adj = np.ceil(5*(N**(0.6)))
    if adj>N+2:
        return N+2
    else:
        return adj

#############################################################################
#############################################################################

def dir2nsim (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        r = re.findall ("_[0-9]+$", dir_name) 
        if ( r ):
            l.append ( int (r[0][1:] ) )
    l.sort() 
    return l

############################################################################
#############################################################################

def dir2temp (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        r = re.findall ("DOP_[0-9]+$", dir_name) 
        if ( r ):
            l.append ( int (r[0][4:] ) )
    l.sort() 
    return l

############################################################################
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

def dir2U (list_of_dirs):
    l = [] 
    for dir_name in list_of_dirs:
        if (re.match("U\d+", dir_name)):
            l.append(dir_name)
    
    l.sort()  
    l = sorted(l, key=lambda x: int(x[1:]) )
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

def get_Rh(master_dict, xlen, ylen, zlen):
    N = master_dict[0][0].shape[0] 
    
    inverse_distance = np.zeros(int(N*(N-1)/2)) 
    tot_steps = len(master_dict) 
    count = 0
    for key in master_dict: 
        coord_arr = unfuck_polymer(master_dict[key][0], xlen, ylen, zlen)
        # print(coord_arr)
        k = 0
        count += 1 
        for i in range(N-1):
            for j in range(i+1,N):
                inverse_distance[k] += 1/(np.linalg.norm( coord_arr[i] - coord_arr[j], 2 ) ) 
                print(coord_arr[i], end=', ')
                print(coord_arr[j], end=', ')
                print( np.linalg.norm( coord_arr[i] - coord_arr[j], 2 ) )
                print ("inverse distance is: " + str(inverse_distance[k]) )
                k += 1
    print ( "sum of inverse distances: " + str (np.sum(inverse_distance) ) )
    print ( "count is " + str(count) )
    ens_average = np.sum(inverse_distance)/count  
    hydrodynamic_radius = 1/ ( ens_average/N**2 )  

    return hydrodynamic_radius 

#############################################################################
#############################################################################

def get_pdict(filename, starting_step, x, y, z):
    f = open(filename, 'r')
    coord_file = f.readlines() 
    
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
    
    for line in coord_file:
        if ( re.search(st_b_str, line) ):
            
            step_num = int ( ( extract_loc_from_string ( line.replace('.', ' ') ) ) )
            
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
            master_dict[step_num][pmer_flag] = np.empty ( (0,3) ) 
            continue 

        elif ( re.search(end_str_1, line) ) and starting_bool: 
            end_step_flag = 1
            step_flag     = 0
            pmer_flag     = -1 
            continue 

        elif ( re.search(end_str_2, line) ) and starting_bool:
            continue 

        elif ( starting_bool ):
            # print ( "step_num is " + str(step_num) )
            monomer_coords                   = extract_loc_from_string ( line ) 
            master_dict[step_num][pmer_flag] = np.vstack ( (master_dict[step_num][pmer_flag], monomer_coords[0:-1] ) ) 
            continue
    
    return master_dict

#############################################################################
#############################################################################

def infiltrate_coords_get_rg ( U, T, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num) 
    print (filename)
    edge = edge_length (dop)
    master_dict = get_pdict (filename, starting_index, edge, edge, edge)

    rg = np.asarray([]) 

    for key in master_dict: 
        rg = np.hstack ( (rg, get_Rg(master_dict[key][0], edge, edge, edge) ) )

    return rg 



#############################################################################
#############################################################################


def infiltrate_coords_get_rh ( U, T, num, dop, coords_files, starting_index ):

    filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/"+ coords_files + "_" + str(num) 
    print (filename)
    edge = edge_length (dop)
    master_dict = get_pdict (filename, starting_index, edge, edge, edge)

    rg = np.asarray([]) 

    for key in master_dict: 
        rg = np.hstack ( (rg, get_Rh(master_dict[key][0], edge, edge, edge) ) )

    return rg 


#############################################################################
#############################################################################


def plot_rg ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
    # dop          = args.dop
        
    if excl_vol_bool:
        U_list.append("Uexcl")

    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = []  
    for U in U_list:
        rg_mean = np.array([]) 
        rg_std  = np.array([])      
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        for T in temperatures:
            
            num_list = np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) )

            pool = multiprocessing.Pool (processes=len(num_list)) 

            print("coords_files = " + coords_files)

            with pool as p:
                results = p.starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), \
                    itertools.repeat(T), num_list, itertools.repeat(dop), itertools.repeat(coords_files), \
                    itertools.repeat(starting_index) ) )
            print(len(results))
            # close the pool 
            pool.close() 
            pool.join() 
            
            rg = np.array([]) 
            for r in results:
                rg = np.hstack ( (rg, r) )
            rg_mean = np.hstack ( ( rg_mean, np.mean(rg) ) )         
            rg_std  = np.hstack ( ( rg_std , np.std(rg)/np.sqrt( len(num_list) ) ) ) 
            
        
        if (U == "Uexcl"):
            ax.errorbar   ( temperatures, rg_mean, yerr=rg_std, linewidth=3 )
        else:
            ax.errorbar   ( temperatures, rg_mean, yerr=rg_std, fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1 
    
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$R_g ^2$", fontsize=18)     
    ax.set_xticks ( np.arange(0, np.max(Tmax)+1, 1 ) )
    plt.savefig   ( "DOP_"+str(dop)+"_rg.png", dpi=1200)
    if show_plot_bool:
        plt.show() 


#############################################################################
#############################################################################


def plot_rg_multi ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
    # dop          = args.dop
        
    if excl_vol_bool:
        U_list.append("Uexcl")

    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = []  
    for U in U_list:
        rg_mean = np.array([]) 
        rg_std  = np.array([])      
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rg_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            num_list = np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) )
            master_num_list.extend (list( num_list ) )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rg_dict[T] = np.asarray([]) 


        # start multiprocessing... 

        pool = multiprocessing.Pool ( processes=len(master_temp_list) ) # len(num_list)) 
        with pool as p:
            results = p.starmap ( infiltrate_coords_get_rg, zip( itertools.repeat(U), master_temp_list, master_num_list, itertools.repeat(dop), itertools.repeat(coords_files), itertools.repeat(starting_index) ) )

        # print (len(results))
        pool.close()
        pool.join() 

         

        for i in range(len(master_temp_list)):
            rg_dict[master_temp_list[i]] = np.hstack ( (rg_dict[master_temp_list[i]], results[i]) )
        
        rg_mean = np.asarray([])
        rg_std  = np.asarray([])  
        for T in temperatures:
            rg_mean = np.hstack ( (rg_mean, np.mean ( rg_dict[T] ) ) )
            rg_std  = np.hstack ( (rg_std , np.std  ( rg_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) )

        if (U == "Uexcl"):
            ax.errorbar   ( temperatures, rg_mean, yerr=rg_std, linewidth=3 )
        else:
            ax.errorbar   ( temperatures, rg_mean, yerr=rg_std, fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1 
    
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$R_g ^2$", fontsize=18)     
    ax.set_xticks ( np.arange(0, np.max(Tmax)+1, 1 ) )
    plt.savefig   ( "DOP_"+str(dop)+"_rg.png", dpi=1200)
    if show_plot_bool:
        plt.show() 



#############################################################################
#############################################################################


def plot_rh ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
    # dop          = args.dop
        
    if excl_vol_bool:
        U_list.append("Uexcl")

    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = []  
    for U in U_list:
        rh_mean = np.array([]) 
        rh_std  = np.array([])      
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        for T in temperatures:
            
            num_list = np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) )

            pool = multiprocessing.Pool (processes=len(num_list)) 

            print("coords_files = " + coords_files)

            with pool as p:
                results = p.starmap ( infiltrate_coords_get_rh, zip( itertools.repeat(U), \
                    itertools.repeat(T), num_list, itertools.repeat(dop), itertools.repeat(coords_files), \
                    itertools.repeat(starting_index) ) )
            print(len(results))
            
            # close the pool 
            pool.close() 
            pool.join() 
            
            rh = np.array([]) 
            for r in results:
                rh = np.hstack ( (rh, r) )
            rh_mean = np.hstack ( ( rh_mean, np.mean(rh) ) )         
            rh_std  = np.hstack ( ( rh_std , np.std(rh)/np.sqrt( len(num_list) ) ) ) 
            
        
        if (U == "Uexcl"):
            ax.errorbar   ( temperatures, rh_mean, yerr=rh_std, linewidth=3, elinewidth=1, capsize=0 )
        else:
            ax.errorbar   ( temperatures, rh_mean, yerr=rh_std, fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1 
    
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$R_h$", fontsize=18)     
    ax.set_xticks ( np.arange(0, np.max(Tmax)+1, 1 ) )
    plt.savefig   ( "DOP_"+str(dop)+"_rg.png", dpi=1200)
    if show_plot_bool:
        plt.show() 


#############################################################################
#############################################################################


def plot_rh_multi ( dop, starting_index, excl_vol_bool, coords_files, show_plot_bool ):

    U_list = dir2U ( os.listdir (".") )
        
    if excl_vol_bool:
        U_list.append("Uexcl")

    # print (edge_length(dop))
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = []  
    for U in U_list:
        rh_mean = np.array([]) 
        rh_std  = np.array([])      
        temperatures = dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        
        # get num_list for each temperature 
        master_temp_list = [] 
        master_num_list = [] 
        rh_dict    = {}
        ntraj_dict = {}
        for T in temperatures: 
            num_list = np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) )
            master_num_list.extend (list( num_list ) )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            rh_dict[T] = np.asarray([]) 


        # start multiprocessing... 

        pool = multiprocessing.Pool ( processes=len(master_temp_list) ) # len(num_list)) 
        with pool as p:
            results = p.starmap ( infiltrate_coords_get_rh, zip( itertools.repeat(U), master_temp_list, \
                master_num_list, itertools.repeat(dop), itertools.repeat(coords_files), \
                itertools.repeat(starting_index) ) )

        # print (len(results))
        pool.close()
        pool.join() 

         

        for i in range(len(master_temp_list)):
            rh_dict[master_temp_list[i]] = np.hstack ( (rh_dict[master_temp_list[i]], results[i]) )
        
        rh_mean = np.asarray([])
        rh_std  = np.asarray([])  
        for T in temperatures:
            rh_mean = np.hstack ( (rh_mean, np.mean ( rh_dict[T] ) ) )
            rh_std  = np.hstack ( (rh_std , np.std  ( rh_dict[T] )/ np.sqrt( ntraj_dict[T] ) ) )

        if (U == "Uexcl"):
            ax.errorbar   ( temperatures, rh_mean, yerr=rh_std, linewidth=3 )
        else:
            ax.errorbar   ( temperatures, rh_mean, yerr=rh_std, fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.copper(i/9), label='_nolegend_' ) 
        i+=1 
    
    ax.legend     ( ["Athermal solvent"], loc='best', fontsize=12)
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$R_g ^2$", fontsize=18)     
    ax.set_xticks ( np.arange(0, np.max(Tmax)+1, 1 ) )
    plt.savefig   ( "DOP_"+str(dop)+"_rg.png", dpi=1200)
    if show_plot_bool:
        plt.show() 


#############################################################################
#############################################################################

def get_flory(U, T, dop_list):
    
    rg_mean = [] 

    for dop in dop_list:
        filename = U+"/DOP_"+str(dop)+"/"+T+"/coords.txt" 
        master_dict = get_pdict (filename, dop+2, dop+2, dop+2)

        rg    = np.asarray([]) 
        steps = np.asarray([]) 
        
        for key in master_dict: 
            rg = np.hstack ( (rg, get_Rg (master_dict[key][0], dop+2, dop+2, dop+2) ) ) 
            steps = np.hstack ( steps, key )
        rg_mean = np.hstack ( rg_mean, np.mean(rg[args.s:] ) )

    rg_mean = np.array (rg_mean)
    
    L = np.polyfit ( np.log(dop_list), np.log(rg_mean), 1, w = np.sqrt(dop_list) ) 
    
    return L[0]/2

#############################################################################
#############################################################################


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
