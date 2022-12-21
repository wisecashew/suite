#!/usr/licensed/anaconda3/2020.7/bin/python

import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import argparse 
import aux
import os
import multiprocessing 
import itertools
import aux 

parser = argparse.ArgumentParser (description="Read a dump file and plot the rdfs.")
parser.add_argument('-T', metavar='T', dest='T', action='store', nargs='+', type=float, help='Enter temperatures you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
args = parser.parse_args() 

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = np.asarray([np.asarray(u) for u in v])

def mysorter_f (x):
	new_str = ""
	for m in x:
		if m.isdigit():
			new_str = new_str+m
	return float(new_str)


if __name__=="__main__":

	start = time.time()
	dop = args.dop
	U_list = aux.dir2U ( os.listdir("."))
	# U_list = ["U1", "U5"]
	fig, ax = plt.subplots(figsize=(8,6))
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort()
	nproc = args.nproc 
	print ("initiate parallel threads...", flush=True)
	pool = multiprocessing.Pool ( processes=nproc ) 
	PLOT_DICT = {}


	for T in temperatures:
		frac_list = [] 
		count_s2s2    = [] 
		master_U_list   = []
		master_num_list = []
		count_s1s1_dict = {}
		count_s2s2_dict = {}
		ntraj_dict = {}

		for U in U_list:
			frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) 
			master_num_list.extend ( num_list )
			master_U_list.extend ( [U]*len (num_list) )
			ntraj_dict [U] = len(num_list)
			count_s2s2_dict [U] = []

		print ("frac_list = ", frac_list, flush=True)
		print ("master_num_list = ",master_num_list, flush=True)
		idx_range = len (master_num_list)//nproc + 1
		print ("idx_range = ", idx_range, flush=True)
		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				print ("in idx = {}".format(u_idx), flush=True)
				print ("about to calculate the s2s2 rdf...", flush=True)
				results_s2s2 = pool.starmap ( aux.extract_rdf_s2s2, zip(master_U_list[u_idx*nproc:], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:]) )
				print ("calculated s2s2 rdf!", flush=True) 
				print ("completed idx = {}".format(u_idx), flush=True) 
			else:
				print ("in idx = {}".format(u_idx), flush=True)
				print ("about to calculate thes2s2 rdf...", flush=True)
				results_s2s2 = pool.starmap (aux.extract_rdf_s2s2, zip(master_U_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc]))
				print ("calculated s2s2 rdf!", flush=True) 
				print ("completed idx = {}".format(u_idx), flush=True) 

			for k in range ( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				count_s2s2_dict[master_U_list[u_idx*nproc+k]].append (results_s2s2[k])

		sorted_U_list = list(np.unique(master_U_list))
		sorted_U_list.sort (key=mysorter_f)
		for U in sorted_U_list:
			count_s2s2.append ( np.mean ( count_s2s2_dict[U] ) )

		PLOT_DICT[T] = np.asarray(count_s2s2)

	pool.close()
	pool.join()

	for T in temperatures:
		print(temperatures)
		f = open ("s2s2.dat", 'w')
		for i in range(len(PLOT_DICT[T])):
			f.write (str(frac_list[i]) + "	" + str(PLOT_DICT[T][i]) + "\n")
		f.close() 

	print ("time = {:.2f}".format(time.time()-start))


