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

def location (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = (lattice_index % (z*z)) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x
	return np.array([xc, yc, zc])

start = time.time()

edge  = 34
dop   = 32
step  = 50000000
frac  = aux.get_frac ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1/U2/DOP_32/0.01/geom_and_esurf.txt")
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

f = open ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1_dump/U2/DOP_32/0.01/lattice_dump_1.mc", 'r')

for line in f:
	if re.findall(start_str, line):
		r = re.findall ("\d+", line)
		step = r[0]
		step_bool = True
		m1_dict[step] = np.zeros((dop,   3))
		s1_dict[step] = np.zeros((nsol1, 3))
		s2_dict[step] = np.zeros((nsol2, 3))
		m_num  = 0
		s1_num = 0
		s2_num = 0
	elif re.findall (end_str, line):
		continue
	elif step_bool:
		info = line.strip().split()
		if info[1] == "m1,":
			m1_dict[step][m_num]  = location(int(info[2]), edge, edge, edge)
			m_num += 1
		elif info[1] == "s1,":
			s1_dict[step][s1_num] = location(int(info[2]), edge, edge, edge)
			s1_num += 1
		elif info[1] == "s2,":
			s2_dict[step][s2_num] = location(int(info[2]), edge, edge, edge)
			s2_num += 1

f.close() 
print ("The dictionaries have been created!", flush=True)
print ("time = {}".format(time.time()-start), flush=True)
print ("time to get pairwise distances...", flush=True)
# i now have all the positions of monomers, solvents and cosolvents
# begin making periodic copies of the solvents 


v = [[0,0,0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, 1, 0], [1, -1, 0], \
[-1, -1, 0], [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1], \
[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], \
[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], \
[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = [np.asarray(u) for u in v]
perturbed_dlist = []
keycount        = 0
keymax          = 150
distances       = np.zeros (dop*nsol1)

delta  = 0.1
rmin   = 0.06
rmax   = edge*np.sqrt(3)
rrange = np.arange (rmin, rmax+2*delta, delta)
hgram  = np.zeros  (len(rrange)-1)
print ("delta = {}".format(delta), flush=True)
print ("rmin  = {}".format(rmin),  flush=True)
print ("rmax  = {}".format(rmax),  flush=True)

start = time.time()
neighbor_count = np.zeros (keymax)

chunk_size = 1000

for key in s1_dict.keys():
	for i in range (27):
		perturbed_dlist.append (s1_dict[key] + edge*v[i])
		print ("periodic image idx = ", i, flush=True)
		for j in range(nsol1//chunk_size + 1):
			start = time.time()
			distances = scipy.spatial.distance_matrix(s1_dict[key], perturbed_dlist[i][(j-1)*(chunk_size):j*chunk_size]).flatten()
			neighbor_count[keycount] += len(distances[ distances < np.sqrt(3)+0.1 ])
			print ("chunk idx = ", j, flush=True)
			end = time.time()
			print ("time for distance_matrix = {:.2f}".format (end-start), flush=True)

		# hist_indices = (distances-rmin)//delta
		# idx, counts  = np.unique(hist_indices, return_counts=True)
		# hgram[idx.astype(int)] += counts
	keycount += 1
	print ("key = {}".format(key), flush=True)
	perturbed_dlist.clear()
	if keycount == keymax:
		break

print ("obtained all distances!") 
print ("time = {}".format(time.time()-start), flush=True)
print ("time to plot...") 

# normalize number of neighbors with number of s2 particles 
neighbor_count = neighbor_count/(nsol2*keymax)

start = time.time()

fig, ax = plt.subplots()
# ax.plot(rrange_clean[:plot_idx+1], hgram_clean[:plot_idx+1], marker='o', markersize=2, color='steelblue')
# ax.axhline(y=0, c='k', linestyle='-', linewidth=0.2)


fig.savefig("rdf.png", dpi=800)

print ("done!")
print ("time = {}".format(time.time()-start), flush=True)

