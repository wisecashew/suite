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
from matplotlib.colors import Normalize
import multiprocessing 
import itertools 
import sys 
import copy
import time
from numba import jit
import scipy.spatial.distance as ssd
from sklearn.linear_model import LinearRegression 

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

def get_info (topology):

    f = open (topology, 'r')
    frac = "frac"
    x    = "x"
    y    = "y"
    z    = "z"
    T    = "kT"
    for line in f:
        if re.findall (frac, line):
            frac_ = re.findall("[0-9]+\.[0-9]+", line)
            continue
        elif re.findall (x, line):
            x_ = re.findall("[0-9]+", line)
            continue
        elif re.findall (y, line):
            y_ = re.findall("[0-9]+", line)
            continue
        elif re.findall (z, line):
            z_ = re.findall("[0-9]+", line)
            continue
        elif re.findall (T, line):
            T_ = re.findall ("[0-9]+\.[0-9]+|[0-9]+\.|[0-9]+", line)
            continue

    f.close()

    return np.array ( [float(x_[0]), float (y_[0]), float (z_[0]), float (T_[0]), float (frac_[0])] )


def get_frac (topology):

    f = open (topology, 'r')
    frac = "frac"
    for line in f:
        if re.findall (frac, line):
            frac_ = re.findall("[0-9]+\.[0-9]+", line)
            break

    f.close()

    return frac_[0]

def get_energy_form_1 (topology):
	f = open (topology, 'r')
	Emm   = "Emm"
	Ems   = "Ems"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"

	for line in f:
		if re.findall (Emm, line):
			r = re.findall (num_re, line)
			# print(r)
			mm = float ( r[0] )
		elif re.findall ( Ems, line):
			r = re.findall (num_re, line)
			# print(r)
			ms = float( r[0] ) 

	return np.array([mm, ms])



def get_energy_form_2 (topology):
	f = open (topology, 'r')
	Emm_a   = "Emm_a"
	Emm_n   = "Emm_n"
	Ems   = "Ems"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"

	for line in f:
		if re.findall (Emm_a, line):
			r = re.findall (num_re, line)
			# print(r)
			mm_a = float ( r[0] )
		elif re.findall (Emm_n, line):
			r = re.findall (num_re, line)
			mm_n = float ( r[0] )
		elif re.findall ( Ems, line):
			r = re.findall (num_re, line)
			# print(r)
			ms = float( r[0] ) 

	return np.array([mm_a, mm_n, ms])

def get_energy_form_3 (topology):
	f = open (topology, 'r')
	Emm_1   = "Emm_1"
	Emm_2   = "Emm_2"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"

	for line in f:
		if re.findall (Emm_1, line):
			r = re.findall (num_re, line)
			mm_1 = float ( r[0] )
		elif re.findall (Emm_2, line):
			r = re.findall (num_re, line)
			mm_2 = float ( r[0] )

	return np.array([mm_1, mm_2])


def get_energy_form_4 (topology):
	f = open (topology, 'r')
	Emm_1   = "Emm_a"
	Emm_2   = "Emm_n"
	B       = "B"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"

	for line in f:
		if re.findall (Emm_1, line):
			r = re.findall (num_re, line)
			mm_1 = float ( r[0] )
		elif re.findall (Emm_2, line):
			r = re.findall (num_re, line)
			mm_2 = float ( r[0] )
		elif re.findall (B, line):
			r = re.findall (num_re, line)
			b = float ( r[0] )

	return np.array([mm_1, mm_2, b])




def get_energy_target (topology):
	f = open (topology, 'r')
	Emm_a   = "Emm_a"
	Emm_n   = "Emm_n"
	Ems1_a  = "Ems1_a"
	Ems1_n  = "Ems1_n"
	Ems2_a  = "Ems2_a"
	Ems2_n  = "Ems2_n"
	Es1s2_a = "Es1s2_a"
	Es1s2_n = "Es1s2_n"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"

	for line in f:
		if re.findall (Emm_a, line):
			r = re.findall (num_re, line)
			# print(r)
			mm_a = float ( r[0] )
		elif re.findall ( Emm_n, line):
			r = re.findall (num_re, line)
			# print(r)
			mm_n = float( r[0] ) 
		elif re.findall ( Ems1_a, line):
			r = re.findall (num_re, line)
			# print(r)
			ms1_a = float( r[0] ) 
		elif re.findall ( Ems1_n, line):
			r = re.findall (num_re, line)
			# print(r)
			ms1_n = float( r[0] ) 
		elif re.findall ( Ems2_a, line):
			r = re.findall (num_re, line)
			# print(r)
			ms2_a = float( r[0] ) 
		elif re.findall ( Ems2_n, line):
			r = re.findall (num_re, line)
			# print(r)
			ms2_n = float( r[0] ) 
		elif re.findall ( Es1s2_a, line):
			r = re.findall (num_re, line)
			# print(r)
			s1s2_a = float( r[0] ) 
		elif re.findall ( Es1s2_n, line):
			r = re.findall (num_re, line)
			# print(r)
			s1s2_n = float(r[0] ) 

	return np.array([mm_a, mm_n, ms1_a, ms1_n, ms2_a, ms2_n, s1s2_a, s1s2_n])


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

def dir2enthalpies (list_of_dirs):
	l = []
	for dir_name in list_of_dirs:
		try:
			l.append (float(dir_name[2:]))
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
        if (re.match("^U\d+$|^U\d+\.\d+$", dir_name)):
            l.append(dir_name)
    
    l.sort()  
    l = sorted(l, key=lambda x: float(x[1:]) )
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

def get_pdict (filename, starting_step, dop, x, y, z):
    
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
            try:
                master_dict[step_num][pmer_flag][m_index] = monomer_coords[0:-1] 
            except KeyError: 
                print ("filename =", filename)
                print ("coords = ",monomer_coords)
                print ("step num =", step_num)
                exit()
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

    master_dict = get_pdict (filename, starting_index, dop, edge, edge, edge)
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

def single_sim_flory_exp_energetic_variation ( U, H, num, dop, coords_file, starting_index, delta ):
	filename = U+"/DOP_"+str(dop)+"/E_"+str(H)+"/"+coords_file+"_"+str(num)+".mc"
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

def get_bond_corr ( U, T, num, dop, coords_file, starting_index, delta):
	filename = U+"/DOP_"+str(dop)+"/"+str(T)+"/"+coords_file+"_"+str(num)+".mc"
	edge     = edge_length  (dop)
	master_dict = get_pdict ( filename, starting_index, dop, edge, edge, edge )
	cos_theta_list = []
	for key in master_dict:
		coord_arr = unfuck_polymer ( master_dict[key][0], edge, edge, edge )
		bonds     = (coord_arr[delta:] - coord_arr [:-delta])
		cos_theta = bonds[1:]*bonds[:-1]/(np.reshape(np.linalg.norm(bonds[1:], axis=1)*np.linalg.norm(bonds[:-1], axis=1), (-1,1) ) )
		cos_theta_list.append ( np.mean (cos_theta) )

	return np.mean (cos_theta_list)


def get_f_of_s ( U, T, num, dop, coords_file, starting_index ):
	
	f_of_s = []
	for i in np.arange (2, 7):
		f_of_s.append ( get_bond_corr ( U, T, num, dop, coords_file, starting_index, int(i) ) )
	
	print (f_of_s)
	ave_l = (1+np.sqrt(2)+np.sqrt(3))/3
	model = LinearRegression()
	model.fit(np.arange(2,7).reshape((-1,1)), np.log(f_of_s))
	r2 = model.score (np.arange(2, 7).reshape((-1,1)), f_of_s)
	
	return (-ave_l/model.coef_, r2)


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
    monomer_order  = np.array([0,0,0])
    magnetization  = []

    for line in f:
        if re.match ( start_str, line ):
            a = re.search ("\d+", line)
            extract_orr = True
            monomer_order = np.array([0,0,0])
            if int ( a.group(0) ) == starting_index:
                start_bool = True
            count          = 0

        elif re.match ( end_str, line ) and start_bool:
            magnetization.append (np.mean (monomer_order))
            monomer_order = np.array ([0,0,0])
            extract_orr = False

        elif extract_orr and start_bool:
            monomer_or = extract_loc_from_string ( line ) [0]
            count = 0
            monomer_order += Or2Dir[monomer_or]

    f.close()

    return np.mean (magnetization)


##########################################################################
##########################################################################



def dump_magnetization ( ortn_file_name, idx, starting_index ):
    Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
            3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
            6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
            9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
            12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
            15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
            18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }
    start_str = "Dumping coordinates at step"
    end_str   = "END"

    f = open(ortn_file_name+"/coords_"+str(idx)+".mc", 'r')

    extract_orr = False
    start_bool  = False
    monomer_order  = np.array([0,0,0])
    magnetization  = []

    for line in f:
        if re.match ( start_str, line ):
            a = re.search ("\d+", line)
            # extract_orr = True
            monomer_order = np.array([0,0,0], dtype=np.float64)
            if int ( a.group(0) ) == starting_index:
                start_bool = True
            count          = 0

        elif re.match ("START", line):
            extract_orr = True

        elif re.match ( end_str, line ) and start_bool:
            magnetization.append (np.linalg.norm(monomer_order))
            monomer_order = np.array ([0,0,0])
            extract_orr = False

        elif extract_orr and start_bool:
            # print (line, )
            # print (extract_loc_from_string(line))
            monomer_or = extract_loc_from_string ( line ) [3]
            count = 0
            monomer_order += Or2Dir[monomer_or]

    f.close()
    step = np.arange (len(magnetization))
    d = {"step": step, "magnetization": magnetization}
    df = pd.DataFrame.from_dict(d)
    df.to_csv (ortn_file_name+"/magnetization_dump_1.mc", index=False, sep='|')

    return



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
        monomer_alignment.clear()
        monomer_monomer   = 0
        # identify all neighbors
        for i in range(N-1):
            for j in range(i+1, N):
                diff = coord_arr[j][:-1]-coord_arr[i][:-1]
                # print ("diff = ",diff)
                if (np.abs(int(diff[0])) == 1 or int(diff[0]) == 0) and (np.abs(int(diff[1])) == 1 or int(diff[1]) == 0) and (np.abs(int(diff[2])) == 1 or int(diff[2]) == 0):
                    monomer_alignment.append( np.dot (Or2Dir[orientations[i]], Or2Dir[orientations[j]]) )
        if len(monomer_alignment) == 0:
            continue
        else:
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

########################################################################################
########################################################################################


def lat2loc (lattice_index, x, y, z):

	R  = np.zeros((len(lattice_index), 3))
	zc = (lattice_index // (z*z))
	yc = ((lattice_index % (z*z)) // y)
	xc = (( ( lattice_index % (z*z) ) % y ) % x)
	R[:,0] = xc
	R[:,1] = yc
	R[:,2] = zc

	return R

def lat2loc_singular (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = ( (lattice_index % (z*z) ) ) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x 

	return np.array ([xc, yc, zc])

def loc2lat (location, x, y, z):
	lat_vec = (location[:,0]%x)+(location[:,1]%y)*y+(location[:,2]%z)*(z*z)
	return lat_vec



########################################################################################
########################################################################################

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = np.asarray([np.asarray(u) for u in v])

'''
def extract_rdf_s1s1 (U, N, T, num):
	dop       = N
	edge      = edge_length (N)
	frac      = get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-dop)*frac))
	nsol1     = edge**3 - dop - nsol2
	rho       = nsol1/edge**3 
	start_str = "FINAL STEP"
	end_str   = "END."
	step_bool = False

	m1_dict = {}
	s1_dict = {}
	s2_dict = {}
	if frac == 1.0:
		return 0
	else: 
		f = open (U + "/DOP_" + str(N) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')

		for line in f:
			if re.findall (start_str, line):
				r = re.findall("\d+", line)
				step = r[0]
				step_bool = True
				m1_dict[step] = np.zeros(dop)
				s1_dict[step] = np.zeros(nsol1)
				s2_dict[step] = np.zeros(nsol2)
				m_num  = 0
				s1_num = 0
				s2_num = 0
			elif re.findall (end_str, line):
				step_bool = False
				m1_dict[step] = np.sort(m1_dict[step], kind='mergesort')
				s1_dict[step] = np.sort(s1_dict[step], kind='mergesort')
				s2_dict[step] = np.sort(s2_dict[step], kind='mergesort') 
				continue 
			elif step_bool:
				info = line.strip().split()
				if info[1] == "m1,":
					m1_dict[step][m_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					m_num += 1
				elif info[1] == "s1,":
					s1_dict[step][s1_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s1_num += 1
				elif info[1] == "s2,":
					s2_dict[step][s2_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s2_num += 1 

		f.close()
		keycount   = 0
		neighbor_count = 0
		solvation_shell = {} 
		for key in s1_dict.keys():
			solvation_shell[key] = np.zeros((0,2))
			for s1_pos in s1_dict[key]:
				neighbors = loc2lat(lat2loc(s1_pos[0], edge, edge, edge) + v, edge, edge, edge).flatten()
				neighbors = neighbors[neighbors <= np.max(s1_dict[key][:,0]]
				x = np.searchsorted(s1_dict[key][:,0], neighbors)
				y = np.equal (s1_dict[key][:,0][x], neighbors)
				solvation_shell[key] = np.vstack((solvation_shell[key], s1_dict[key][y]))
			u, indices = np.unique(solvation_shell[key][:,0], return_indices=True)
			solvation_shell[key] = solvation_shell[key][indices]
			keycount += 1
		final_count = neighbor_count/(nsol1*keycount)

	return final_count 
'''
########################################################################################
########################################################################################

'''
def extract_rdf_s2s2 (U, N, T, num):
	dop       = N
	edge      = edge_length (N)
	frac      = get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-dop)*frac))
	nsol1     = edge**3 - dop - nsol2
	rho       = nsol1/edge**3 
	start_str = "FINAL STEP"
	end_str   = "END."
	step_bool = False

	m1_dict = {}
	s1_dict = {}
	s2_dict = {}
	if frac == 0.0:
		return 0
	else:
		f = open (U + "/DOP_" + str(N) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')
		for line in f:
			if re.findall (start_str, line):
				r = re.findall("\d+", line)
				step = r[0]
				step_bool = True 
				m1_dict[step] = np.zeros(dop)
				s1_dict[step] = np.zeros(nsol1)
				s2_dict[step] = np.zeros(nsol2)
				m_num  = 0
				s1_num = 0
				s2_num = 0
			elif re.findall (end_str, line):
				step_bool = False 
				m1_dict[step] = np.sort(m1_dict[step], kind='mergesort')
				s1_dict[step] = np.sort(s1_dict[step], kind='mergesort')
				s2_dict[step] = np.sort(s2_dict[step], kind='mergesort') 
				continue 
			elif step_bool:
				info = line.strip().split()
				if info[1] == "m1,":
					m1_dict[step][m_num] = int(info[2])  # location (int(info[2]), edge, edge, edge)
					m_num += 1
				elif info[1] == "s1,":
					s1_dict[step][s1_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s1_num += 1
				elif info[1] == "s2,":
					s2_dict[step][s2_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s2_num += 1 

		f.close()
		keycount   = 0
		neighbor_count = 0
		for key in s2_dict.keys():
			for s2_pos in s2_dict[key]:
				neighbors = loc2lat(lat2loc(s2_pos, edge, edge, edge)+v, edge, edge, edge).flatten()
				x = np.searchsorted(s2_dict[key], neighbors)
				neighbor_count += np.sum( np.hstack((s2_dict[key], -1))[x]==neighbors )
			keycount += 1
		final_count = neighbor_count/(keycount*nsol2)
	return final_count
'''
###############################################################################################
###############################################################################################
'''
def extract_rdf_s1s2 (U, N, T, num):
	dop       = N
	edge      = edge_length (N)
	frac      = get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-dop)*frac))
	nsol1     = edge**3 - dop - nsol2
	rho       = nsol1/edge**3 
	start_str = "FINAL STEP"
	end_str   = "END."
	step_bool = False

	m1_dict = {}
	s1_dict = {}
	s2_dict = {}
	f = open (U + "/DOP_" + str(N) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')
	for line in f:
		if re.findall (start_str, line):
			r = re.findall("\d+", line)
			step = r[0]
			step_bool = True 
			m1_dict[step] = np.zeros(dop)
			s1_dict[step] = np.zeros(nsol1)
			s2_dict[step] = np.zeros(nsol2)
			m_num  = 0
			s1_num = 0
			s2_num = 0
		elif re.findall (end_str, line):
			step_bool = False 
			m1_dict[step] = np.sort(m1_dict[step], kind='mergesort')
			s1_dict[step] = np.sort(s1_dict[step], kind='mergesort')
			s2_dict[step] = np.sort(s2_dict[step], kind='mergesort') 
			continue 
		elif step_bool:
			info = line.strip().split()
			if info[1] == "m1,":
				m1_dict[step][m_num] = int(info[2])  # location (int(info[2]), edge, edge, edge)
				m_num += 1
			elif info[1] == "s1,":
				s1_dict[step][s1_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
				s1_num += 1
			elif info[1] == "s2,":
				s2_dict[step][s2_num] = int(info[2]) # location (int(info[2]), edge, edge, edge)
				s2_num += 1 

	f.close()
	keycount   = 0
	neighbor_count = 0
	if nsol1 == 0 or nsol2 == 0:
		final_count = 0
	else:
		for key in s1_dict.keys():
			for s1_pos in s1_dict[key]:
				neighbors = loc2lat(lat2loc(s1_pos, edge, edge, edge)+v, edge, edge, edge).flatten()
				x = np.searchsorted(s2_dict[key], neighbors)
				neighbor_count += np.sum( np.hstack((s2_dict[key], -1))[x]==neighbors )
			keycount += 1
			if keycount == 20:
				break
		final_count = neighbor_count/(keycount*nsol1*nsol2/(edge**3))
	# print ("final_count = ",final_count, flush=)
	keycount_   = 0
	neighbor_count_ = 0
	
	if nsol1 == 0 or nsol2 == 0:
		final_count_ = 0
	else:
		for key in s2_dict.keys():
			for s2_pos in s2_dict[key]:
				neighbors = loc2lat(lat2loc(s2_pos, edge, edge, edge)+v, edge, edge, edge).flatten()
				x = np.searchsorted(s1_dict[key], neighbors)
				neighbor_count_ += np.sum( np.hstack((s1_dict[key], -1))[x]==neighbors )
			keycount_ += 1
			if keycount_ == 20:
				break
		final_count_ = neighbor_count_/(keycount_*nsol2)
	
	return final_count # (final_count+final_count_)/2
'''
###############################################################################################
###############################################################################################
def latvec2loc (latvec, edge):
	s = list(np.shape(latvec))
	s.append(3)
	s = tuple(s)
	w = np.zeros(s)
	w[:,2] = (x[:] // (edge*edge)).reshape(s[0],)
	w[:,1] = (x[:] % (edge*edge) // edge).reshape(s[0],)
	w[:,0] = ( ( ( x[:] % (edge*edge) ) % edge ) % edge).reshape(s[0],)
	return w 


Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
            3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
            6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
            9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
            12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
            15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
            18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }


def applyOr2Dir(num):
	return Or2Dir[num]

def s1s2_alignment_check ( s1, s2, edge ):

	r = lat2loc (s2[:,0], edge, edge, edge)-lat2loc_singular (s1[0], edge, edge, edge) # this is the connection vector
	# take dot product of r with orientation of s1 
	t1 = np.arccos( np.clip( ( np.dot( r/np.linalg.norm (r, axis=1)[:,None], Or2Dir[s1[1]] ) ), -1, 1) )
	Or = np.array  ( list( map( applyOr2Dir, s2[:,1]) ) )
	t2 = np.arccos ( np.clip( np.sum (Or*-r, axis=1), -1, 1))
	summation = t1+t2
	aligned = np.sum( summation < np.pi/2)
	misaligned = np.sum ( summation >= np.pi/2)
	return aligned, misaligned 


def extract_s1s2_aligned_solvation_shell (U, N, T, num):

	start = time.time()
	dop       = N
	edge      = edge_length (N)
	frac      = get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-dop)*frac))
	nsol1     = edge**3 - dop - nsol2
	rho       = nsol1/edge**3 
	start_str = "FINAL STEP"
	end_str   = "END."
	step_bool = False
	np.printoptions(Suppress=True)
	m1_dict = {}
	s1_dict = {}
	s2_dict = {}
	
	if frac == 1.0 or frac == 0.0:
		return 0
	else:
		# print ("About to obtain the dictionaries for all particles...", flush=True)
		f = open (U + "/DOP_" + str(N) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')
		for line in f:
			if re.findall (start_str, line):
				r = re.findall("\d+", line)
				step = r[0]
				step_bool = True 
				m1_dict[step] = np.zeros((dop  , 2))
				s1_dict[step] = np.zeros((nsol1, 2))
				s2_dict[step] = np.zeros((nsol2, 2))
				m_num  = 0
				s1_num = 0
				s2_num = 0
			elif re.findall (end_str, line):
				step_bool = False 
				m1_dict[step] = np.sort(m1_dict[step], axis=0, kind='mergesort')
				s1_dict[step] = np.sort(s1_dict[step], axis=0, kind='mergesort')
				s2_dict[step] = np.sort(s2_dict[step], axis=0, kind='mergesort') 
				continue 
			elif step_bool:
				info = line.strip().split()
				if info[1] == "m1,":
					m1_dict[step][m_num][0] = int(info[2])  # location (int(info[2]), edge, edge, edge)
					m1_dict[step][m_num][1] = int(info[0][:-1])
					m_num += 1
				elif info[1] == "s1,":
					s1_dict[step][s1_num][0] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s1_dict[step][s1_num][1] = int(info[0][:-1])
					s1_num += 1
				elif info[1] == "s2,":
					s2_dict[step][s2_num][0] = int(info[2]) # location (int(info[2]), edge, edge, edge)
					s2_dict[step][s2_num][1] = int(info[0][:-1])
					s2_num += 1
		f.close()
	# print ("Obtained dictionaries for all types of particles!", flush=True)
	# obtain all atoms in the solvation shell of type s1 
	
	# print ("Time to find the solvation shell for the polymer at each snapshot...", flush=True)
	solvation_shell = {}
	for key in m1_dict.keys():
		R = lat2loc (m1_dict[key][:,0], edge, edge, edge)
		R = R[:, None] + v
		R = R.reshape(R.shape[0]*R.shape[1], 3) % edge
		neighbors = loc2lat (R, edge, edge, edge).flatten()
		x = np.searchsorted( s1_dict[key][:,0], neighbors)
		neighbors = neighbors [ x < len(s1_dict[key][:,0]) ]
		x = x [ x < len(s1_dict[key][:,0]) ] 
		y = np.equal ( s1_dict[key][:,0][x], neighbors )
		s1_particles = s1_dict[key][x][y]
		u, indices = np.unique (s1_particles[:,0], return_index=True)
		solvation_shell[key] = s1_particles[indices]
	
	# print ("Found the solvation shell for the polymer at each snapshot!", flush=True)
	
	# print ("Counting the number of aligned and misaligned s1-s2 interactions...", flush=True)
	aligned_count = 0
	misaligned_count = 0
	for key in solvation_shell.keys():
		for s1 in solvation_shell[key]:
			neighbors = loc2lat(lat2loc_singular(s1[0], edge, edge, edge)+v, edge, edge, edge).flatten()
			x         = np.searchsorted(s2_dict[key][:,0],neighbors)
			neighbors = neighbors [ x < len(s2_dict[key][:,0]) ]
			x         = x [ x < len(s2_dict[key][:,0]) ]
			y         = np.equal (s2_dict[key][:,0][x], neighbors)
			if np.sum(y) == 0:
				aligned = 0
				misaligned = 0
			else:
				aligned, misaligned = s1s2_alignment_check ( s1, s2_dict[key][x][y], edge)
			aligned_count    += aligned
			misaligned_count += misaligned

	# print ("# of aligned interactions = ",aligned_count, flush=True)
	# print ("Counted the number of aligned and misaligned s1-s2 interactions!", flush=True)
	# print ("time = {:.2f} seconds.".format(time.time()-start), flush=True)
	return aligned_count

######################################################################
######################################################################

class PiecewiseNormalize(Normalize):
    def __init__(self, xvalues, cvalues):
        self.xvalues = xvalues
        self.cvalues = cvalues

        Normalize.__init__(self)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        if self.xvalues is not None:
            x, y = self.xvalues, self.cvalues
            return np.ma.masked_array(np.interp(value, x, y))
        else:
            return Normalize.__call__(self, value, clip)


def gradient_image(ax, extent, direction, cmap_range=(0, 1), **kwargs):
    # divnorm = matplotlib.colors.TwoSlopeNorm (vcenter=0.33, vmin=0.2, vmax=0.8)
    divnorm = PiecewiseNormalize ([0.2, 0.3, 0.7, 0.8], [0.1, 0.2, 0.9, 1.0])
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    extent
        The extent of the image as (xmin, xmax, ymin, ymax).
        By default, this is in Axes coordinates but may be
        changed using the *transform* keyword argument.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular useful is *cmap*.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, extent=extent, interpolation='bicubic', norm=divnorm, **kwargs)
    return im
