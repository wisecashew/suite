#!/usr/licensed/anaconda3/2020.7/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/current/Explicit_Solvation/py_analysis')
import aux 
import pandas as pd
import argparse
import itertools
import multiprocessing
import os

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.stdout.flush()

parser = argparse.ArgumentParser(description="Go into lattice dump and get contacts.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-T', dest='T', action='store', nargs='+', type=float, help='Enter enthalpies you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('--csv-name', type=str, dest='cn', action='store', help='Enter name of database.')
args = parser.parse_args() 

######################################################

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], \
[0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
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
def get_contacts (U, dop, T, num):

	
	df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(T)+"/energydump_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
	mm_contacts    = df["mm_tot"].values[-2000:]
	mma_contacts   = df["mm_aligned"].values[-2000:]
	mmn_contacts   = df["mm_naligned"].values[-2000:]
	ms_contacts    = df["ms1_tot"].values[-2000:] + df["ms2_tot"].values[-2000:]
	ms1_contacts   = df["ms1_tot"].values[-2000:]
	ms1a_contacts  = df["ms1_aligned"].values[-2000:]
	ms1n_contacts  = df["ms1_naligned"].values[-2000:]
	ms2_contacts   = df["ms2_tot"].values[-2000:]
	ms2a_contacts  = df["ms2_aligned"].values[-2000:]
	ms2n_contacts  = df["ms2_naligned"].values[-2000:]
	s1s2_contacts  = df["ms1s2_tot"].values[-2000:]
	s1s2a_contacts = df["ms1s2_aligned"].values[-2000:]
	s1s2n_contacts = df["ms1s2_naligned"].values[-2000:]

	return np.array([np.mean(mm_contacts), np.mean(mma_contacts), np.mean(mmn_contacts), np.mean(ms_contacts), np.mean(ms1_contacts), np.mean(ms1a_contacts), np.mean (ms1n_contacts), np.mean(ms2_contacts), np.mean(ms2a_contacts), np.mean(ms2n_contacts), np.mean(s1s2_contacts), np.mean(s1s2a_contacts), np.mean(s1s2n_contacts)])

#######################################################################
#######################################################################


if __name__=="__main__":
	
	start = time.time()
	U_list = aux.dir2U (os.listdir("."))
	print (U_list)
	dop = args.dop
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort()
	pd.options.display.float_format = '{:,.2f}'.format
	i = 0
	nproc = args.nproc
	pool= multiprocessing.Pool (processes=nproc)
	CONTACTS_DICT = {}

	for T in temperatures:
		CONTACTS_DICT["T"]       = []
		CONTACTS_DICT["x"]       = []
		CONTACTS_DICT["M1-M1"]   = []
		CONTACTS_DICT["M1-M1-A"] = []
		CONTACTS_DICT["M1-M1-N"] = []
		CONTACTS_DICT["M1-S"]    = []
		CONTACTS_DICT["M1-S1"]   = []
		CONTACTS_DICT["M1-S1-A"] = []
		CONTACTS_DICT["M1-S1-N"] = []
		CONTACTS_DICT["M1-S2"]   = []
		CONTACTS_DICT["M1-S2-A"] = []
		CONTACTS_DICT["M1-S2-N"] = []
		CONTACTS_DICT["S1-S2"]   = []
		CONTACTS_DICT["S1-S2-A"] = []
		CONTACTS_DICT["S1-S2-N"] = []
		print ( "We are in T = "+str(T)+", and N = " + str(dop) + "...", flush=True)
		frac_list         = []
		master_U_list     = []
		# master_index_list = []
		master_num_list   = []
		ntraj_dict        = {}
		contacts_dict     = {}
		for U in U_list:
			frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt") )
			num_list = list ( np.unique ( aux.dir2nsim (os.listdir (U + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
			master_num_list.extend ( num_list )
			master_U_list.extend ([U]*len(num_list))
			ntraj_dict[U] = len(num_list)
			contacts_dict [U] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_U_list

		idx_range = len(master_num_list)//nproc + 1
		for u_idx in range(idx_range):
			if u_idx == idx_range-1:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:] ) )
			else:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc] ) )
			
			print ( "Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

			for k in range( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				contacts_dict[master_U_list[u_idx*nproc+k]].append (results[k])

		sorted_U_list = U_list # list (np.unique (master_U_list))
		# sorted_U_list.sort(key=mysorter_f)
		CONTACTS_DICT["T"].extend( [T]*len(sorted_U_list) )
		CONTACTS_DICT["x"].extend( frac_list )
		for U in sorted_U_list:
			CONTACTS_DICT["M1-M1"].append ( np.mean( np.array (contacts_dict[U])[:,0] ) )
			CONTACTS_DICT["M1-M1-A"].append ( np.mean( np.array (contacts_dict[U])[:,1] ) )
			CONTACTS_DICT["M1-M1-N"].append ( np.mean( np.array (contacts_dict[U])[:,2] ) )
			CONTACTS_DICT["M1-S"].append ( np.mean( np.array (contacts_dict[U])[:,3] ) )
			CONTACTS_DICT["M1-S1"].append ( np.mean( np.array (contacts_dict[U])[:,4] ) )
			CONTACTS_DICT["M1-S1-A"].append ( np.mean( np.array (contacts_dict[U])[:,5] ) )
			CONTACTS_DICT["M1-S1-N"].append ( np.mean( np.array (contacts_dict[U])[:,6] ) )
			CONTACTS_DICT["M1-S2"].append ( np.mean( np.array (contacts_dict[U])[:,7] ) )
			CONTACTS_DICT["M1-S2-A"].append ( np.mean( np.array (contacts_dict[U])[:,8] ) )
			CONTACTS_DICT["M1-S2-N"].append ( np.mean( np.array (contacts_dict[U])[:,9] ) )
			CONTACTS_DICT["S1-S2"].append ( np.mean( np.array (contacts_dict[U])[:,10] ) )
			CONTACTS_DICT["S1-S2-A"].append ( np.mean( np.array (contacts_dict[U])[:,11] ) )
			CONTACTS_DICT["S1-S2-N"].append ( np.mean( np.array (contacts_dict[U])[:,12] ) )

	pool.close()
	pool.join ()

	df = pd.DataFrame.from_dict (CONTACTS_DICT)
	df.to_csv (args.cn, sep='|', index=False)
	stop = time.time()
	print ("Time to get the contacts database is " + str(stop-start) + " seconds.", flush=True)



