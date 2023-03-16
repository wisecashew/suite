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

parser = argparse.ArgumentParser (description="Read a dump file and plot the rdfs.")
parser.add_argument('-T', metavar='T', dest='T', action='store', nargs='+', type=float, help='Enter temperatures you want plotted.')
parser.add_argument('-nproc', metavar='N', type=int, dest='nproc', action='store', help='Request these many proccesses.')
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('--png-name', metavar='image name', dest='pn', type=str, action='store', help='enter image name.')
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
	fig, ax = plt.subplots(figsize=(8,6))
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort()
	nproc = args.nproc 
	print ("initiate parallel threads...", flush=True)
	pool = multiprocessing.Pool ( processes=nproc ) 
	PLOT_DICT = {}
	# print (temperatures)

	for T in temperatures:
		frac_list = [] 
		count_s1s2    = []
		master_U_list   = []
		master_num_list = []
		count_s1s2_dict = {}
		ntraj_dict = {}

		for U in U_list:
			frac_list.append ( aux.get_frac(U+"/geom_and_esurf.txt"))
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) 
			master_num_list.extend ( num_list )
			master_U_list.extend ( [U]*len (num_list) )
			ntraj_dict [U] = len(num_list)
			count_s1s2_dict [U] = [] 

		# print ("frac_list = ", frac_list, flush=True)
		# print ("master_num_list = ",master_num_list, flush=True)
		idx_range = len (master_num_list)//nproc + 1
		# print ("idx_range = ", idx_range, flush=True)
		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				print ("in idx = {}".format(u_idx), flush=True)
				results_s1s2 = pool.starmap (aux.extract_rdf_s1s2, zip(master_U_list[u_idx*nproc:], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:]) )
				print ("obtained the s1s2 rdf!", flush=True)
				print ("completed idx = {}".format(u_idx), flush=True) 
			else:
				print ("in idx = {}".format(u_idx), flush=True)
				results_s1s2 = pool.starmap (aux.extract_rdf_s1s2, zip(master_U_list[u_idx*nproc:(u_idx+1)*nproc], itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc]))
				print ("obtained the s1s2 rdf!", flush=True)
				print ("completed idx = {}".format(u_idx), flush=True) 

			for k in range ( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				count_s1s2_dict[master_U_list[u_idx*nproc+k]].append (results_s1s2[k])

		sorted_U_list = list(np.unique(master_U_list))
		sorted_U_list.sort (key=mysorter_f)
		print ("sorted U list = ", sorted_U_list)
		for U in sorted_U_list:
			count_s1s2.append ( np.mean ( count_s1s2_dict[U] ) )
		PLOT_DICT[T] = np.asarray (count_s1s2)

	pool.close()
	pool.join()
	
	i = 0
	if len(temperatures) == 1:
		colors = cm.coolwarm (np.linspace(0,1,3))
		i=1
	else:
		colors = cm.coolwarm (np.linspace(0,1,len(temperatures)))

	for T in temperatures:
		ax.plot (frac_list, PLOT_DICT[T]/26, linestyle='-', marker='o', markeredgecolor='k', linewidth=3, c=colors[i], label='_nolegend', markersize=10)
	i += 1 



	plt.savefig (args.pn, bbox_inches='tight', dpi=1200)
	print ("time = {:.2f}".format(time.time()-start))

