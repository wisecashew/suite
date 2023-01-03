#!/usr/licensed/anaconda3/2020.7/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import aux 
import pandas as pd
import argparse
import itertools
import multiprocessing
import os
import sys

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

parser = argparse.ArgumentParser(description="Go into lattice dump and get contacts.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-T', metavar='T', dest='T', action='store', nargs='+', type=float, help='Enter temperatures you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
args = parser.parse_args() 

######################################################

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = [np.asarray(u) for u in v]

#######################################################

def get_starting_ind ( U, T, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + str(T) + "/" + dumpfile + "_" + str(num) + ".mc"
	df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	L = len(df["energy"])
	return int(df["time_step"].values[L-2000])

def mysorter_f (x):
	new_str = ""
	for m in x:
		if m.isdigit():
			new_str = new_str+m
	return float(new_str)

def lat2loc (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = (lattice_index % (z*z)) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x

	return np.array([int(xc), int(yc), int(zc)])

def loc2lat (location, x, y, z):
	lat_vec = (location[:,0]%x)+(location[:,1]%y)*y+(location[:,2]%z)*(z*z)
	return lat_vec


# create a dictionary of all locations 
def get_contacts (U, dop, T, num, starting_index):

	all_dict = {}
	edge      = aux.edge_length (dop)
	frac      = aux.get_frac(U+"/DOP_"+str(dop)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-dop)*frac))
	nsol1     = edge**3 - dop - nsol2
	rho       = nsol1/edge**3
	start_str = "FINAL STEP:"
	end_str   = "END."
	step_bool = False

	# print ("Obtaining a dictionary for the trajectory information...", flush=True)

	f = open (U + "/DOP_" + str(dop) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')
	for line in f:
		if re.findall (start_str, line):
			r = re.findall("\d+", line)
			step =int( r[0] )
			if step < starting_index:
				continue
			else:
				step_bool = True 
				all_dict[step] = list(np.zeros(dop+nsol1+nsol2))
				all_num = 0
		elif re.findall (end_str, line) and step_bool:
			step_bool = False
			continue
		elif step_bool:
			info = line.strip().split()
			all_dict[step][int(info[2])] = info[1][:-1]
			all_num += 1

	f.close()

	# print ("The master dictionary has been made for the entire trajectory!", flush=True)
	# print ("Obtaining contacts...", flush=True)

	contact_dict = {} 
	contact_dict [("m1", "m1")] = []
	contact_dict [("m1", "s1")] = []
	contact_dict [("m1", "s2")] = []
	contact_dict [("s1", "s1")] = []
	contact_dict [("s1", "s2")] = []
	contact_dict [("s2", "s2")] = []
	contact_dict ["timestep"]   = []

	# dictionary has been found
	# now find all contacts

	i = 0
	for key in all_dict.keys():
		# print ("Step = ", key, flush=True)
		for k in contact_dict:
			if k == "timestep":
				contact_dict[k].append(key)
			else:
				contact_dict[k].append(0)
		lat_idx = 0
		
		for particle in all_dict[key]:
			neighbors = loc2lat (lat2loc (lat_idx, edge, edge, edge) + v, edge, edge, edge).flatten()
			x = [all_dict[key][neigh] for neigh in neighbors]
			m1_count = x.count("m1")
			s1_count = x.count("s1")
			s2_count = x.count("s2")
			t = [tuple(sorted([particle, "m1"])), tuple(sorted([particle, "s1"])), tuple(sorted([particle, "s2"]))]
			contact_dict[t[0]][i] += m1_count/2
			contact_dict[t[1]][i] += s1_count/2
			contact_dict[t[2]][i] += s2_count/2
			lat_idx += 1
		i += 1

	df = pd.DataFrame.from_dict (contact_dict, orient='columns') 
	df.columns = ["m1-m1", "m1-s1", "m1-s2", "s1-s1", "s1-s2", "s2-s2", "timestep"]
	# print (df)
	return np.array([np.mean(df["m1-m1"].values), np.mean(df["m1-s1"].values), np.mean(df["m1-s2"]), np.mean(df["s1-s1"]), np.mean(df["s1-s2"]), np.mean(df["s2-s2"])])

#######################################################################
#######################################################################


if __name__=="__main__":
	
	start = time.time()
	U_list = aux.dir2U (os.listdir("."))
	dop = args.dop
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort()

	i = 0
	nproc = args.nproc
	pool= multiprocessing.Pool (processes=nproc)
	CONTACTS_DICT = {}
	CONTACTS_DICT["T"] = []
	CONTACTS_DICT["x"] = []
	CONTACTS_DICT["M1-M1"] = []
	CONTACTS_DICT["M1-S1"] = []
	CONTACTS_DICT["M1-S2"] = []
	CONTACTS_DICT["S1-S1"] = []
	CONTACTS_DICT["S1-S2"] = []
	CONTACTS_DICT["S2-S2"] = []

	for T in temperatures:
		print ( "We are in T = "+str(T)+", and N = " + str(dop) + "...", flush=True)
		frac_list         = []
		master_U_list     = []
		master_index_list = []
		master_num_list   = []
		ntraj_dict        = {}
		contacts_dict     = {}
		for U in U_list:
			frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt") )
			num_list = list ( np.unique ( aux.dir2nsim (os.listdir (U + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
			master_num_list.extend ( num_list )
			for num in num_list:
				master_index_list.append ( get_starting_ind (U, T, num, dop, "energydump") )
			master_U_list.extend ([U]*len(num_list))
			ntraj_dict[U] = len(num_list)
			contacts_dict [U] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_U_list

		idx_range = len(master_num_list)//nproc + 1
		for u_idx in range(idx_range):
			if u_idx == idx_range-1:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:], master_index_list[u_idx*nproc:] ) )
			else:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc], master_index_list[u_idx*nproc:] ) )
			
			print ( "Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )
		
			for k in range( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				contacts_dict[master_U_list[u_idx*nproc+k]].append (results[k])

		sorted_U_list = list (np.unique (master_U_list))
		sorted_U_list.sort(key=mysorter_f)
		CONTACTS_DICT["T"].extend( [T]*len(sorted_U_list) )
		CONTACTS_DICT["x"].extend( frac_list )
		for U in sorted_U_list:
			CONTACTS_DICT["M1-M1"].append ( np.mean( np.array (contacts_dict[U])[:,0] ) )
			CONTACTS_DICT["M1-S1"].append ( np.mean( np.array (contacts_dict[U])[:,1] ) )
			CONTACTS_DICT["M1-S2"].append ( np.mean( np.array (contacts_dict[U])[:,2] ) )
			CONTACTS_DICT["S1-S1"].append ( np.mean( np.array (contacts_dict[U])[:,3] ) )
			CONTACTS_DICT["S1-S2"].append ( np.mean( np.array (contacts_dict[U])[:,4] ) )
			CONTACTS_DICT["S2-S2"].append ( np.mean( np.array (contacts_dict[U])[:,5] ) )


	pool.close()
	pool.join ()

	df = pd.DataFrame.from_dict (CONTACTS_DICT)
	df.to_csv ("CONTACTS-COSOLVENT_T_"+str(temperatures[0])+".csv", sep='|', index=False)
	stop = time.time()
	print ("Time to get the contacts database is " + str(stop-start) + " seconds.", flush=True)




