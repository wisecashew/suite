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
parser.add_argument('-U', metavar='UX', dest='U', type=str, action='store', help='Enter a potential energy surface.')
parser.add_argument('-T', metavar='T', dest='T', type=float, action='store', help='Enter a temperature.') 
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('--png-name', dest='pn', type=str, action='store', help='enter a name for the png.')
args = parser.parse_args() 

v = [[0,0,0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = [np.asarray(u) for u in v] 

def location (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = (lattice_index % (z*z)) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x
	return np.array([xc, yc, zc])

def extract_rdf (U, N, T, num):
	
	edge      = aux.edge_length (N)
	frac      = aux.get_frac(args.U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
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
			m1_dict[step] = np.zeros((dop,   3))
			s1_dict[step] = np.zeros((nsol1, 3))
			s2_dict[step] = np.zeros((nsol2, 3))
			m_num  = 0
			s1_num = 0
			s2_num = 0
		elif re.findall (end_str, line):
			step_bool = False 
			continue 
		elif step_bool:
			info = line.strip().split()
			if info[1] == "m1,":
				m1_dict[step][m_num] = location (int(info[2]), edge, edge, edge)
				m_num += 1
			elif info[1] == "s1,":
				s1_dict[step][s1_num] = location (int(info[2]), edge, edge, edge)
				s1_num += 1
			elif info[1] == "s2,":
				s2_dict[step][s2_num] = location (int(info[2]), edge, edge, edge)
				s2_num += 1 

	f.close() 
	perturbed_dlist = [] 
	keycount = 0
	distances = np.zeros (dop*nsol1) 
	delta  = 0.05 
	rmin   = 0.01 
	rmax   = 10 # edge*np.sqrt(3)/2
	rrange = np.arange (rmin, rmax+2*delta, delta)
	hgram  = np.zeros  (len(rrange)-1)

	for key in m1_dict.keys():
		for i in range(27):
			perturbed_dlist.append (s1_dict[key] + edge*v[i])
			distances = scipy.spatial.distance_matrix (m1_dict[key], perturbed_dlist[i]).flatten()
			distances = distances[ distances < rmax+delta ]
			hist_indices = (distances-rmin)//delta 
			idx, counts  = np.unique (hist_indices, return_counts=True)
			hgram[idx.astype(int)] += counts 
		keycount += 1
		perturbed_dlist.clear()
		# if keycount == 20:
		# 	break

	hgram = hgram / (dop*keycount) 
	hgram = hgram / (4*np.pi*((rrange[:-1]+delta/2)**2)*delta*rho)

	return (hgram, rrange[:-1])


if __name__=="__main__":
	start = time.time()
	dop = args.dop
	U   = args.U
	T   = args.T 
	num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) 
	# num_list = num_list[0:1]
	nproc = len(num_list)
	print ("initiate parallel threads...", flush=True)
	pool = multiprocessing.Pool ( processes=nproc ) 
	results = pool.starmap (extract_rdf, zip(itertools.repeat(U), itertools.repeat(dop), itertools.repeat(T), num_list))
	print ("Pool closed. Pool had {} threads.".format (len(results)), flush=True)
	pool.close()
	pool.join() 
	
	
	hgram = np.zeros (len(results[0][0]))
	for i in range(nproc):
		hgram += results[i][0] 
	hgram  = hgram/nproc
	rrange = results[0][1]

	for i in range(len(hgram)):
		if hgram[i] != 0.:
			zidx = i
			break 

	# hgram_clean  = hgram[0:zidx]
	# hgram_clean  = np.hstack((hgram_clean, hgram[hgram>0]))
	# rrange_clean = rrange[0:zidx]
	# rrange_clean = np.hstack((rrange_clean, rrange[hgram>0]))
	hgram_clean  = hgram
	rrange_clean = rrange
	# from scipy.signal import savgol_filter 
	# yhat = savgol_filter (hgram_clean[zidx:], 51, 3)
	fig, ax = plt.subplots()
	ax.plot(rrange_clean, hgram_clean, color='steelblue', markersize=2, marker='o')
	# ax.plot(rrange_clean[zidx:], yhat, c='coral')
	# ax.plot(rrange_clean[:zidx], hgram_clean[:zidx], c='coral')
	# ax.plot(rrange_clean, np.hstack((hgram_clean[:zidx],yhat)), c='coral')
	w = 20
	hgram_con = np.convolve(hgram_clean, np.ones(w), 'valid')/w
	rrange_con = np.convolve(rrange_clean, np.ones(w), 'valid')/w
	ax.plot(rrange_con, hgram_con, color='coral', markersize=2, marker='o')
	ax.axhline (y=1, c='k', linestyle='-')
	ax.set_ylim(bottom=-0.5)
	fig.savefig (args.pn, dpi=1200)
	stop = time.time()
	print ("time = {} seconds.".format(stop-start))

