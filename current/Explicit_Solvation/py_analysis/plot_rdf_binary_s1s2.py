#!~/.conda/envs/data_analysis/bin/python


import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
import aux 

def lat2loc (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = (lattice_index % (z*z)) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x

	return np.array([int(xc), int(yc), int(zc)])

def loc2lat (location, x, y, z):
	lat_vec = (location[:,0]%x)+(location[:,1]%y)*y+(location[:,2]%z)*(z*z)
	return lat_vec

v = [ [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1] ]
v = np.asarray([np.asarray(u) for u in v])

start = time.time()
U = "U5"
edge  = 34
dop   = 32
step  = 50000000
frac  = aux.get_frac ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1_dump/"+U+"/DOP_32/0.01/geom_and_esurf.txt")

nsol2 = int(np.floor(((edge**3)-dop)*frac))
nsol1 = edge**3 - dop - nsol2 
rho   = nsol1/edge**3
start_str = "FINAL STEP" 
end_str = "END."

step_bool = False

m1_dict = {}
s1_dict = {}
s2_dict = {}

print ("initializations complete. infiltrate dump files!", flush=True)

f = open ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1_dump/"+U+"/DOP_32/0.01/lattice_dump_1.mc", 'r')

for line in f:
	if re.findall(start_str, line):
		r = re.findall ("\d+", line)
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
			m1_dict[step][m_num]  = int(info[2]) # location(int(info[2]), edge, edge, edge)
			m_num += 1
		elif info[1] == "s1,":
			s1_dict[step][s1_num] = int(info[2]) # location(int(info[2]), edge, edge, edge)
			s1_num += 1
		elif info[1] == "s2,":
			s2_dict[step][s2_num] = int(info[2]) # location(int(info[2]), edge, edge, edge)
			s2_num += 1

f.close() 
print ("The dictionaries have been created!", flush=True)
print ("time = {}".format(time.time()-start), flush=True)
print ("time to get pairwise distances...", flush=True)
print ("number of solvent one molecules = ", nsol1, flush=True)
print ("number of solvent two molecules = ", nsol2, flush=True)
print ("frac = ",frac, flush=True)
# i now have all the positions of monomers, solvents and cosolvents
# begin making periodic copies of the solvents 

keycount        = 0
keymax          = 50
neighbor_count  = 0
start = time.time()

for key in s1_dict.keys():
	for s1_pos in s1_dict[key]:
		neighbors = loc2lat(lat2loc(s1_pos, edge, edge, edge) + v, edge, edge, edge).flatten()
		x = np.searchsorted(s1_dict[key], neighbors)
		neighbor_count += np.sum( np.hstack((s1_dict[key], -1))[x]==neighbors )
	keycount += 1
	print ("key = {}".format(key), flush=True)
	if keycount == keymax:
		break

neighbor_count = neighbor_count/(nsol1*keycount)
print ("neighbor_count = ", neighbor_count)
print ("obtained all distances!") 
print ("time = {}".format(time.time()-start), flush=True)
print ("time to plot...") 

fig, ax = plt.subplots()


print ("done!")
print ("time = {}".format(time.time()-start), flush=True)

