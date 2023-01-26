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
parser.add_argument('--png-name', dest='pn', type=str, action='store', help='enter a name for the png.')
args = parser.parse_args() 

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = [np.asarray(u) for u in v]

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


def extract_rdf_s1s1 (U, N, T, num):
	
	edge      = aux.edge_length (N)
	frac      = aux.get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
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
		for key in s1_dict.keys():
			for s1_pos in s1_dict[key]:
				neighbors = loc2lat(lat2loc(s1_pos, edge, edge, edge) + v, edge, edge, edge).flatten()
				x = np.searchsorted(s1_dict[key], neighbors)
				neighbor_count += np.sum( np.hstack((s1_dict[key], -1))[x] == neighbors )
			keycount += 1
		final_count = neighbor_count/(nsol1*keycount)

	return final_count 


def extract_rdf_s2s2 (U, N, T, num):

	edge      = aux.edge_length (N)
	frac      = aux.get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
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

if __name__=="__main__":

	start = time.time()
	dop = args.dop
	U_list = aux.dir2U ( os.listdir("."))
	U_list = ["U5"]
	fig, ax = plt.subplots(figsize=(8,6))
	temperatures = [float(elem) for elem in args.T]
	temperatures.sort()
	nproc = args.nproc 
	print ("initiate parallel threads...", flush=True)
	pool = multiprocessing.Pool ( processes=nproc ) 
	PLOT_DICT = {}


	for T in temperatures:
		frac_list = [] 
		count_s1s1    = []
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
			count_s1s1_dict [U] = [] 
			count_s2s2_dict [U] = []

		print ("frac_list = ", frac_list, flush=True)
		print ("master_num_list = ",master_num_list, flush=True)
		idx_range = len (master_num_list)//nproc + 1
		print ("idx_range = ", idx_range, flush=True)
		for u_idx in range (idx_range):
			if u_idx == idx_range-1:
				print ("in idx = {}".format(u_idx), flush=True)
				results_s1s1 = pool.starmap (extract_rdf_s1s1, zip(itertools.repeat(U), itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:]) )
				print ("obtained the s1s1 rdf!", flush=True)
				print ("about to calculate thes2s2 rdf...", flush=True)
				results_s2s2 = pool.starmap (extract_rdf_s2s2, zip(itertools.repeat(U), itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:]) ) 
				print ("completed idx = {}".format(u_idx), flush=True) 
			else:
				print ("in idx = {}".format(u_idx), flush=True)
				results_s1s1 = pool.starmap (extract_rdf_s1s1, zip(itertools.repeat(U), itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc]))
				print ("obtained the s1s1 rdf!", flush=True)
				print ("about to calculate thes2s2 rdf...", flush=True)
				results_s2s2 = pool.starmap (extract_rdf_s2s2, zip(itertools.repeat(U), itertools.repeat(dop), itertools.repeat(T), master_num_list[u_idx*nproc:(u_idx+1)*nproc]))
				print ("completed idx = {}".format(u_idx), flush=True) 

			for k in range ( len( master_U_list[u_idx*nproc:(u_idx+1)*nproc] ) ):
				count_s1s1_dict[master_U_list[u_idx*nproc+k]].append (results_s1s1[k])
				count_s2s2_dict[master_U_list[u_idx*nproc+k]].append (results_s2s2[k])

		sorted_U_list = list(np.unique(master_U_list))
		sorted_U_list.sort (key=mysorter_f)
		for U in sorted_U_list:
			count_s1s1.append ( np.mean ( count_s1s1_dict[U] ) )
			count_s2s2.append ( np.mean ( count_s2s2_dict[U] ) )


		PLOT_DICT[T] = (np.asarray (count_s1s1), np.asarray(count_s2s2) ) 

	pool.close()
	pool.join()

	for T in temperatures:
		print(temperatures)
		f = open ("s1s1.dat", 'w')
		f.write ("x	count	\n")
		for i in range(len(PLOT_DICT[T][0])):
			f.write (str(frac_list[i]) + "	" + str(PLOT_DICT[T][0][i]) + "\n")
		f.close()

		f = open ("s2s2.dat", 'w')
		for i in range(len(PLOT_DICT[T][1])):
			f.write (str(frac_list[i]) + "	" + str(PLOT_DICT[T][1][i]) + "\n")
		f.close() 

		ax.bar ( frac_list, PLOT_DICT[T][0], color='crimson', width=0.05)
		ax.bar ( frac_list, PLOT_DICT[T][1], bottom=PLOT_DICT[T][0], color='darkorchid', width=0.05)

	ax.axhline (y=0, c='k', linewidth=0.2)
	plt.savefig (args.pn, dpi=1000)
	print ("time = {:.2f}".format(time.time()-start))


