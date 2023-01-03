#!~/.conda/envs/data_analysis/bin/python


import time
import numpy as np
import re 
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')
import scipy 
import scipy.spatial
# import aux 

def location (lattice_index, x, y, z):
	zc = lattice_index // (z*z)
	yc = (lattice_index % (z*z)) // y
	xc = ( ( lattice_index % (z*z) ) % y ) % x
	return np.array([xc, yc, zc])

start = time.time()

edge  = 34
dop   = 32
step  = 50000000
frac  = 0 # aux.get_frac ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1/U1/DOP_32/0.01/geom_and_esurf.txt")
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

f = open ("/scratch/gpfs/satyend/MC_POLYMER/SingleChainComprehensive/runs/all_regimes/cosolvent/isotropic1_dump/U1/DOP_32/0.01/lattice_dump_1.mc", 'r')

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

for key in m1_dict.keys():
	for i in range (27):
		perturbed_dlist.append (s1_dict[key] + edge*v[i])
		distances = scipy.spatial.distance_matrix(m1_dict[key], perturbed_dlist[i]).flatten()
		distances = distances[ distances < rmax+delta ]
		hist_indices = (distances-rmin)//delta
		idx, counts  = np.unique(hist_indices, return_counts=True)
		hgram[idx.astype(int)] += counts
	keycount += 1
	print ("key = {}".format(key), flush=True)
	perturbed_dlist.clear()
	if keycount == 20:
		break

print ("obtained all distances!") 
print ("time = {}".format(time.time()-start), flush=True)
print ("time to plot...") 

start = time.time()
hgram = hgram / (dop*keycount) 
hgram = hgram / (4*np.pi*((rrange[:-1]+delta/2)**2)*delta*rho) 

# hunt for the first nonzero index 
for i in range(len(hgram)):
	if hgram[i] != 0.:
		iidx = i
		break 

hgram_clean  = hgram[0:iidx]
hgram_clean  = np.hstack((hgram_clean, hgram[hgram>0]))
rrange_clean = rrange[0:iidx]
rrange_clean = np.hstack((rrange_clean, rrange[:-1][hgram>0]))

plot_idx = np.where(rrange_clean>0.9*edge)[0][0]

fig, ax = plt.subplots()
ax.plot(rrange_clean[:plot_idx+1], hgram_clean[:plot_idx+1], marker='o', markersize=2, color='steelblue')
ax.axhline(y=1, c='k', linestyle='-')

def smooth (y, box_pts):
	box = np.ones(box_pts)/box_pts 
	y_smooth = np.convolve (y, box, mode='same')
	return y_smooth 

filter_r = rrange[:-1][hgram>0]
filter_h = hgram[hgram>0]
print (filter_h[0:10])

ax.plot (filter_r[:len(filter_r)//2-100], smooth(filter_h[:len(filter_h)//2-100], 2), c='coral')

from scipy.signal import savgol_filter 
yhat = savgol_filter (filter_h[:len(filter_h)//2-100], 51, 3)
ax.plot (filter_r[:len(filter_r)//2-100], yhat, c='pink')
ax.plot (np.hstack((rrange[0:iidx],1)), np.hstack((hgram[0:iidx],yhat[0])), c='pink')
print (hgram[iidx])
# yhat   = savgol_filter (hgram_clean[:len(hgram_clean)//2-100], 51, 3)
# ax.plot (rrange_clean[:len(rrange_clean)//2-100], yhat, c='pink')

fig.savefig("rdf.png", dpi=800)

print ("done!")
print ("time = {}".format(time.time()-start), flush=True)

