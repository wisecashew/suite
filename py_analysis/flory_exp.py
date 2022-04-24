#!/usr/licensed/anaconda3/2020.7/bin/python
import numpy as np
import re 
import os 
import matplotlib.pyplot as plt

st = "stat_kt_"

dop = [8, 16, 24, 32, 48, 64]
txtfiles = []

for d in dop:
	txtfiles.append(st+str(d)+".txt") 

rg = [] 
rg_std = [] 
for file in txtfiles: 
	f = open(file, 'r') 
	data = f.readlines()
	f.close() 
	m = data[1].strip().split() 
	rg.append(float(m[-1]))
	q = data[3].strip().split() 
	rg_std.append( float(q[-1]) )

plt.errorbar(x=np.log10(dop), y=np.log10(rg), yerr=rg_std, fmt='-.', ecolor='lightgreen', elinewidth=1, capsize=5 )
model = np.polyfit(np.log10(dop), np.log10(rg), 1) 
pred = model[0]*np.log10(dop)+model[1]
plt.plot(np.log10(dop), pred, 'r', alpha=0.3)
plt.xlabel("log($N_D$)")
plt.ylabel("log$\langle R_g \\rangle$")
plt.title("Scaling of radius of gyration with degree of polymerization") 
string = "$\\nu$ = {:0.2f}".format(model[0])
plt.text(1.4,0.3,string)
plt.savefig("flory_exp_T_0.5.png", dpi=1200)
plt.show()

j = open('rg_with_dop.txt', 'w') 
for l in range(len(dop)): 
	j.write(str(dop[l]) + " | " + str(rg[l]) +"\n" ) 
j.close()
