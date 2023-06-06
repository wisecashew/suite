#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description = "Read a LAMMPS .avg file for gyration.")
parser.add_argument('-i', metavar='f1.avg f2.avg ...', dest='i', nargs='+', help='enter address of xvg file')
parser.add_argument('-j', metavar='f1.avg f2.avg ...', dest='j', nargs='+', help='enter address of xvg file')
parser.add_argument('-T', nargs='+', dest='T', metavar='T1 T2 ..', help='temperature of each simulation') 
args = parser.parse_args()

rg1 = np.asarray([])

fig = plt.figure( figsize=(1.75,1.4) )
plt.margins(x=0, y=0)
ax = plt.axes()
ax.tick_params (direction='in', bottom=True, top=True, left=True, right=True, pad=2)
ax.set_ylim (bottom=0, top=30)
ax.set_xlim (0, 140)
# ax.minorticks_on()
ax.set_yticks(np.arange(0,31, 5))
ax.set_xticks(np.arange(0,160,20))
# ax.xaxis.set_tick_params(which='minor', direction='in', right=True)
# ax.yaxis.set_tick_params(which='minor', direction='in', right=True)

# ax.xaxis.grid (True, which='minor')

# ax.yaxis.grid (True, which='minor')
# plt.minorticks_on()
plt.style.use('classic')
for fi in args.i:
    with open(fi, 'r') as f:
        for line in f:
            if line.startswith('# TimeStep c_rog[0] c_rog[1] c_rog[2] c_rog[3]'): 
                break
        df = pd.read_csv(f, names=["t", "Rg", "Rgx", "Rgy", "Rgz"], delim_whitespace=True)
        rg1 = np.hstack( (rg1, df["Rg"]) )
        f.close()

plt.plot( 5000*(np.arange(20, len(rg1)+20))[:27000:100]/10**6 , rg1[:27000:100], color=(0.8, 0.1, 0.1), linewidth=1)
# plt.xlabel("t (ns)", fontsize=4, labelpad=2)
# plt.ylabel("$R_g (\\AA)$", fontsize=4, labelpad=2)

rg2 = np.asarray([])

if args.j == None:
    pass
else:
    for fi in args.j:
        with open(fi, 'r') as f:
            for line in f:
                if line.startswith('# TimeStep c_rog[0] c_rog[1]'):
                    break
            df = pd.read_csv(f, names=["t", "Rg", "Rgx", "Rgy", "Rgz"], delim_whitespace=True)
            rg2 = np.hstack ( (rg2, df["Rg"]) ) 
            f.close()

plt.plot( 5000*np.arange(20,len(rg2)+20)[:27000:100]/10**6, rg2[:27000:100], linewidth=1, color=(0.1, 0.1, 0.8) )
plt.xticks (fontsize=5)
plt.yticks (fontsize=5)

# plt.legend(args.T)
# print( "Mean radius of gyration is {:0.2f}".format( np.mean(rg1) ))  
plt.tight_layout()
plt.savefig( "rg_multiplot.png", dpi=1200 )
plt.show()
