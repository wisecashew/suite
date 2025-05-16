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
parser.add_argument('--U', dest='U', action='store', nargs='+', type=str, help='Enter enthalpies you want plotted.')
parser.add_argument('--frac', dest='frac', action='store', nargs='+', type=float, help='Enter fracs you want probed.')
parser.add_argument('-s', dest='s', action='store', type=int, help='Enter when to start collecting contacts.')
parser.add_argument('--dump', dest='dump', action='store', type=str, help='Enter name of contacts dump file.')
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

def get_starting_ind ( U, f, num, dop, dumpfile):
	filename = U + "/DOP_" + str(dop) + "/" + str(f) + "/" + dumpfile + "_" + str(num) + ".mc"
	df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(f)+"/"+args.dump+"_"+str(num)+".mc", sep=' \| ', names=["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "timestep"], engine='python', skiprows=1)
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
def get_contacts (U, dop, f, num, s):

	df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(f)+"/"+args.dump+"_"+str(num)+".mc", sep=' \| ', names=["total", "total_s", "total_c", "aligned_sc", "misaligned_sc", "timestep"], engine='python', skiprows=1)
	total         = df["total"].values[-s:]
	total_s       = df["total_s"].values[-s:]
	total_c       = df["total_c"].values[-s:]
	aligned_sc    = df["aligned_sc"].values[-s:]
	misaligned_sc = df["misaligned_sc"].values[-s:]

	return np.array([total, total_s, total_c, aligned_sc, misaligned_sc])

#######################################################################
#######################################################################


if __name__=="__main__":
	
	start     = time.time()
	U_list    = args.U
	dop       = args.dop
	frac_list = args.frac
	s         = args.s
	i         = 0
	nproc     = args.nproc
	pool      = multiprocessing.Pool (processes=nproc)
	CONTACTS_DICT      = {}
	CONTACTS_DICT["U"] = []
	CONTACTS_DICT["x"] = []
	CONTACTS_DICT["total"]   = []
	CONTACTS_DICT["total_s"] = []
	CONTACTS_DICT["total_c"] = []
	CONTACTS_DICT["aligned_sc"]    = []
	CONTACTS_DICT["misaligned_sc"]   = []

	for U in U_list:
		print ( "We are in U = "+str(U)+", and N = " + str(dop) + "...", flush=True)
		master_frac_list     = []
		master_num_list   = []
		ntraj_dict        = {}
		contacts_dict     = {}
		for f in frac_list:
			num_list = list (np.unique (aux.dir2nsim (os.listdir (U + "/DOP_" + str(dop) + "/" + str(f) ), args.dump ) ) )
			master_frac_list.extend([f]*len(num_list))
			master_num_list.extend ( num_list )
			ntraj_dict    [f] = len(num_list)
			contacts_dict [f] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_U_list

		idx_range = len(master_num_list)//nproc + 1
		for u_idx in range(idx_range):
			if u_idx == idx_range-1:
				results = pool.starmap ( get_contacts, zip( itertools.repeat(U), itertools.repeat(dop), master_frac_list[u_idx*nproc:], master_num_list[u_idx*nproc:], itertools.repeat(s) ) )
			else:
				results = pool.starmap ( get_contacts, zip( itertools.repeat(U), itertools.repeat(dop), master_frac_list[u_idx*nproc:(u_idx+1)*nproc], master_num_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(s) ) )
			
			print ( "Pool has been closed. This pool had {} threads.".format (len(results) ), flush=True )
		
			for k in range( len( master_frac_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				contacts_dict[master_frac_list[u_idx*nproc+k]].append (results[k])

		CONTACTS_DICT["U"].extend([U]*len(frac_list))
		CONTACTS_DICT["x"].extend( frac_list )
		# print(contacts_dict)
		for f in frac_list:
			CONTACTS_DICT["total"].append   ( np.mean( np.array (contacts_dict[f])[:,0] ) )
			CONTACTS_DICT["total_s"].append ( np.mean( np.array (contacts_dict[f])[:,1] ) )
			CONTACTS_DICT["total_c"].append ( np.mean( np.array (contacts_dict[f])[:,2] ) )
			CONTACTS_DICT["aligned_sc"].append    ( np.mean( np.array (contacts_dict[f])[:,3] ) )
			CONTACTS_DICT["misaligned_sc"].append   ( np.mean( np.array (contacts_dict[f])[:,4] ) )

	pool.close()
	pool.join ()

	df = pd.DataFrame.from_dict (CONTACTS_DICT)
	df.to_csv (args.cn, sep='|', index=False)
	stop = time.time()
	print ("Time to get the contacts database is " + str(stop-start) + " seconds.", flush=True)

