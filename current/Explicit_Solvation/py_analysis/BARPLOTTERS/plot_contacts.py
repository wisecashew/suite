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


# create a dictionary of all locations 
def get_contacts():
	all_dict = {}
	contact_dict = {} 
	U = "U1"
	N = 32
	T = "1.0"
	dop = 32
	num = 1
	edge      = aux.edge_length (N)
	frac      = aux.get_frac(U+"/DOP_"+str(N)+"/"+str(T)+"/geom_and_esurf.txt")
	nsol2     = int(np.floor(((edge**3)-N)*frac))
	nsol1     = edge**3 - N - nsol2
	rho       = nsol1/edge**3 
	starting_index = 25000000
	start_str = "FINAL STEP:"
	end_str   = "END."
	step_bool = False

	print ("Obtaining a dictionary for the trajectory information...", flush=True)

	f = open (U + "/DOP_" + str(N) + "/" + str(T) + "/lattice_dump_" + str(num) + ".mc", 'r')
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
			# all_dict[step] = np.sort  (all_dict[step], kind='mergesort')
			continue
		elif step_bool:
			info = line.strip().split()
			# print (info)
			all_dict[step][int(info[2])] = info[1][:-1]
			all_num += 1

	f.close()

	print ("The master dictionary has been made for the entire trajectory!", flush=True)
	print ("Obtaining contacts...", flush=True)

	contact_dict = {} 
	contact_dict [("m1", "m1")] = []
	contact_dict [("m1", "s1")] = []
	contact_dict [("m1", "s2")] = []
	contact_dict [("s1", "s1")] = []
	contact_dict [("s1", "s2")] = []
	contact_dict [("s2", "s2")] = []
	contact_dict ["timestep"] = []
	# dictionary has been found
	# now find all contacts
	i = 0
	for key in all_dict.keys():
		print ("Step = ", key, flush=True)
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
	df.to_csv ("COSOLVENT-CONTACTS.csv", sep='|', index=False)
	print ("average, m1-m1:", np.mean(df["m1-m1"].values),", m1-s1:", np.mean(df["m1-s1"].values), ", m1-s2:", np.mean(df["m1-s2"].values), ", s1-s1:", np.mean(df["s1-s1"].values), ", s1-s2: ",np.mean(df["s1-s2"].values), ", s2-s2:", np.mean(df["s2-s2"].values) )
	print ("contacts = ", np.mean(df["m1-m1"].values)+np.mean(df["m1-s1"].values)+np.mean(df["m1-s2"].values)+np.mean(df["s1-s1"].values)+np.mean(df["s1-s2"].values)+np.mean(df["s2-s2"].values))

if __name__=="__main__":
	print ("Execute...", flush=True)
	start = time.time()
	get_contacts()
	stop = time.time()
	print ("Time = ", stop-start)

