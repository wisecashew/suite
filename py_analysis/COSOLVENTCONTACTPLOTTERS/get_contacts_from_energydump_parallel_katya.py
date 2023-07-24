#!/home/satyend/.conda/envs/phase/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
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
parser.add_argument('-s', dest='s', action='store', type=int, help='Enter when to start collecting contacts.')
parser.add_argument('--csv-name', dest='cn', action='store', type=str, help='Enter name of CSV file to create.')
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

def get_starting_ind ( U, H, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + H + "/" + dumpfile + "_" + str(num) + ".mc"
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
def get_contacts (U, frac, num, s):

	df = pd.read_csv(str(U)+"/"+frac+"/energydump_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	mm_contacts    = df["mm_tot"].values[-s:]
	mma_contacts   = df["mm_aligned"].values[-s:]
	mmn_contacts   = df["mm_naligned"].values[-s:]
	ms_contacts    = df["ms1_tot"].values[-s:]+df["ms2_tot"].values[-s:]
	ms1_contacts   = df["ms1_tot"].values[-s:]
	ms1a_contacts  = df["ms1_aligned"].values[-s:]
	ms1n_contacts  = df["ms1_naligned"].values[-s:]
	ms2_contacts   = df["ms2_tot"].values[-s:]
	ms2a_contacts  = df["ms2_aligned"].values[-s:]
	ms2n_contacts  = df["ms2_naligned"].values[-s:]
	s1s2_contacts  = df["ms1s2_tot"].values[-s:]
	s1s2a_contacts = df["ms1s2_aligned"].values[-s:]
	s1s2n_contacts = df["ms1s2_naligned"].values[-s:]

	return np.array([np.mean(mm_contacts), np.mean(mma_contacts), np.mean(mmn_contacts), np.mean(ms_contacts), np.mean(ms1_contacts), np.mean(ms1a_contacts), np.mean(ms1n_contacts), np.mean(ms2_contacts), np.mean(ms2a_contacts), np.mean(ms2n_contacts), np.mean(s1s2_contacts), np.mean(s1s2a_contacts), np.mean(s1s2n_contacts)])

#######################################################################
#######################################################################


if __name__=="__main__":
	
	start = time.time()
	U_list = aux.dir2u (os.listdir("."))
	dop = args.dop
	s = args.s
	i = 0
	nproc = args.nproc
	pool= multiprocessing.Pool (processes=nproc)
	CONTACTS_DICT            = {}
	CONTACTS_DICT["u"]       = []
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

	for u in U_list:
		print ( "We are in U = " + str(u) + "...", flush=True )
		frac_list         = list (os.listdir (u))
		xu_list           = [u]*len(frac_list)
		frac_list.sort ()
		master_U_list     = []
		master_frac_list  = []
		master_num_list   = []
		ntraj_dict        = {}
		contacts_dict     = {}
		contacts_dict.clear ()

		for frac in frac_list:
			num_list = list (np.unique (aux.dir2nsim (os.listdir (u + "/" + frac) ) ) )
			master_num_list.extend  ( num_list )
			master_frac_list.extend ( [frac]*len(num_list) )
			master_U_list.extend    ( [u]*len(num_list) )
			ntraj_dict    [frac] = len ( num_list )
			contacts_dict [frac] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_U_list

		idx_range = len(master_num_list) // nproc + 1
		for u_idx in range(idx_range):
			if u_idx == idx_range-1:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:], master_frac_list[u_idx*nproc:], master_num_list[u_idx*nproc:], itertools.repeat(s) ) )
			else:
				results = pool.starmap ( get_contacts, zip( master_U_list[u_idx*nproc:(u_idx+1)*nproc], master_frac_list[u_idx*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(s) ) )
			
			print ( "Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )

			for k in range( len( master_frac_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				contacts_dict[master_frac_list[u_idx*nproc+k]].append (results[k])

		CONTACTS_DICT["u"].extend(xu_list)
		# CONTACTS_DICT["x"].extend(frac_list)
		for frac in frac_list:
			CONTACTS_DICT["x"].append (float(frac[-3:]))

		for frac in frac_list:
			CONTACTS_DICT["M1-M1"].append   ( np.mean( np.array (contacts_dict[frac])[:,0]  ) )
			CONTACTS_DICT["M1-M1-A"].append ( np.mean( np.array (contacts_dict[frac])[:,1]  ) )
			CONTACTS_DICT["M1-M1-N"].append ( np.mean( np.array (contacts_dict[frac])[:,2]  ) )
			CONTACTS_DICT["M1-S"].append    ( np.mean( np.array (contacts_dict[frac])[:,3]  ) )
			CONTACTS_DICT["M1-S1"].append   ( np.mean( np.array (contacts_dict[frac])[:,4]  ) )
			CONTACTS_DICT["M1-S1-A"].append ( np.mean( np.array (contacts_dict[frac])[:,5]  ) )
			CONTACTS_DICT["M1-S1-N"].append ( np.mean( np.array (contacts_dict[frac])[:,6]  ) )
			CONTACTS_DICT["M1-S2"].append   ( np.mean( np.array (contacts_dict[frac])[:,7]  ) )
			CONTACTS_DICT["M1-S2-A"].append ( np.mean( np.array (contacts_dict[frac])[:,8]  ) )
			CONTACTS_DICT["M1-S2-N"].append ( np.mean( np.array (contacts_dict[frac])[:,9]  ) )
			CONTACTS_DICT["S1-S2"].append   ( np.mean( np.array (contacts_dict[frac])[:,10] ) )
			CONTACTS_DICT["S1-S2-A"].append ( np.mean( np.array (contacts_dict[frac])[:,11] ) )
			CONTACTS_DICT["S1-S2-N"].append ( np.mean( np.array (contacts_dict[frac])[:,12] ) )

	pool.close()
	pool.join ()
	# print (CONTACTS_DICT)

	df = pd.DataFrame.from_dict (CONTACTS_DICT)
	df.to_csv (args.cn, sep='|', index=False)
	stop = time.time()
	print ("Time to get the contacts database is " + str(stop-start) + " seconds.", flush=True)


