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
parser.add_argument('-dop',       metavar='DOP', dest='dop',     type=int,     action='store', help='enter a degree of polymerization.')
parser.add_argument('--U',        dest='U',      action='store', nargs='+',    type=str,       help='Enter enthalpies you want plotted.')
parser.add_argument('--T',        dest='T',      action='store', nargs='+',    type=float,     help='Enter temperatures you want probed.')
parser.add_argument('--csv-name', dest='cn',     action='store', type=str,     help='Enter name of CSV file to create.')
parser.add_argument('-nproc',     metavar='N',   type=int,       dest='nproc', action='store', help='Request these many proccesses.')
args = parser.parse_args() 

######################################################

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = np.array(v) # np.asarray(u) for u in v]

#######################################################

def get_starting_ind ( U, f, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + str(f) + "/" + dumpfile + "_" + str(num) + ".mc"
	df = pd.read_csv(filename, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0)
	L = len(df["energy"])
	return 0 #int(df["time_step"].values[L-2000])

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
def get_contacts (U, dop, f, num, s):

	df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(f)+"/energydump_"+str(num)+".mc", sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"], engine='python', skiprows=0) 
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

def count_solvent(pdict):
	count = []
	for idx in pdict:
		# print(f"@idx = {idx}:", flush=True, end=' ')
		polymer     = pdict[idx][0]
		result      = np.array(polymer[:, np.newaxis, :] + v[np.newaxis, :, :], dtype=int)
		result      = result.reshape(-1, 3)
		result      = np.unique(result, axis=0)
		set_polymer = {tuple(row) for row in polymer}
		set_result  = {tuple(row) for row in result}
		diff_set    = set_result - set_polymer
		result      = np.array(list(result))
		# print(f"{len(result)}", flush=True)
		count.append(len(result))
	return count
	
#######################################################################

def get_solvation_shell (U, dop, f, num):

	filename = U + "/DOP_" + str(dop) + "/" + str(f) + "/coords_" + str(num)+".mc" 
	# print(f"filename = {filename}")
	edge            = aux.edge_length (dop)
	starting_index  = get_starting_ind(U, f, num, dop, "energydump")
	master_dict     = aux.get_pdict (filename, starting_index, dop, edge, edge, edge)
	# print(master_dict)
	solvation_shell = count_solvent(master_dict)
	return solvation_shell

#######################################################################

if __name__=="__main__":
	
	start     = time.time()
	U_list    = args.U
	dop       = args.dop
	T_list    = args.T
	i         = 0
	nproc     = args.nproc
	pool      = multiprocessing.Pool (processes=nproc)
	CONTACTS_DICT              = dict()
	CONTACTS_DICT["U"]         = []
	CONTACTS_DICT["T"]         = []
	CONTACTS_DICT["timestep"]  = []
	CONTACTS_DICT["SOLVATION"] = []

	# print(get_solvation_shell("F10", 32, "0.0", 1, s))
	# exit()

	for U in U_list:
		print ( "We are in U = "+str(U)+", and N = " + str(dop) + "...", flush=True)
		master_T_list     = []
		master_num_list   = []
		ntraj_dict        = {}
		contacts_dict     = {}
		for t in T_list:
			num_list = ["1"] # list (np.unique (aux.dir2nsim (os.listdir (U + "/DOP_" + str(dop) + "/" + str(t) ) ) ) )
			master_T_list.extend([t]*len(num_list))
			master_num_list.extend ( num_list )
			ntraj_dict    [t] = len(num_list)
			contacts_dict [t] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_U_list
		idx_range = len(master_num_list)//nproc + 1
		for u_idx in range(idx_range):
			if u_idx == idx_range-1:
				solv    = pool.starmap ( get_solvation_shell, zip( itertools.repeat(U), itertools.repeat(dop), master_T_list[u_idx*nproc:], master_num_list[u_idx*nproc:] ) )
			else:
				solv = pool.starmap ( get_solvation_shell, zip( itertools.repeat(U), itertools.repeat(dop), master_T_list[u_idx*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc]) )
			
			print("Pool has been closed. This pool had {} threads.".format(len(solv)), flush=True)
			for k in range(len(master_T_list[u_idx*nproc:(u_idx+1)*nproc])):
				contacts_dict[master_T_list[u_idx*nproc+k]].append (solv[k])

		for t in T_list:
			print(len(contacts_dict[t]))
			CONTACTS_DICT["U"].extend([U]*len(contacts_dict[t][0]))
			CONTACTS_DICT["T"].extend([t]*len(contacts_dict[t][0]))
			CONTACTS_DICT["SOLVATION"].extend(contacts_dict[t][0])
			CONTACTS_DICT["timestep"].extend(list(np.arange(0,len(CONTACTS_DICT["SOLVATION"]))))

	pool.close()
	pool.join ()

	df = pd.DataFrame.from_dict (CONTACTS_DICT)
	df.to_csv (args.cn, sep='|', index=False)
	stop = time.time()
	print ("Time to get the contacts database is " + str(stop-start) + " seconds.", flush=True)

