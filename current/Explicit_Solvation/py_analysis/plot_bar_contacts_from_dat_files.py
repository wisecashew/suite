#!/usr/licensed/anaconda3/2020.7/bin/python

import re
import pandas as pd 
import argparse 
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib
from matplotlib.ticker import StrMethodFormatter
matplotlib.use('Agg') 
import matplotlib.cm as cm
import matplotlib.ticker as tck
import argparse 
import aux 
import os 

frac_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
df = pd.read_csv("m1sn.dat", header=0, sep='\t')
cols = df.keys()
fig, ax = plt.subplots(figsize=(8,6))

ax.bar (frac_list, df[cols[0]].values, color='coral', width=0.05, edgecolor='k')
ax.bar (frac_list, df[cols[1]].values, bottom=df[cols[0]].values, color='steelblue', width=0.05, edgecolor='k')
ax.bar (frac_list, df[cols[2]].values, bottom=df[cols[0]].values+df[cols[1]].values, color='seagreen', width=0.05, edgecolor='k')

ax.legend(["$m_1$-$s_1$ contacts", "$m_1$-$s_2$ contacts", "$m_1$-$m_1$ contacts"], ncol=3)
ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
ax.tick_params ( axis='x', labelsize=16, direction="in", left="off", labelleft="on" )
ax.tick_params ( axis='y', labelsize=16, direction="in", left="off", labelleft="on" )
# cbar = plt.colorbar(sm, orientation='vertical')
# cbar.set_ticks ( [0, 0.25, 0.5, 0.75, 1.0] )
# cbar.set_ticklabels( ["$10^{-2}$", "$10^{-1}$", "1", "$10^{1}$", "$10^2$"] ) 
# cbar.set_ticks( [] )
# cbar.ax.tick_params (labelsize=14)
# cbar.set_ticklabels( [] )
# cbar.ax.set_ylabel ( "Strength of better solvent \n", fontsize=16, rotation=270 ) 
# ax.set_xscale('log')
# ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) ) 
ax.axhline (y=0, c='k', linewidth=0.2)
# ax.yaxis.get_major_locator().set_params(integer=True)
ax.minorticks_on()
ax.set_xlim((-0.05, 1.05))
ax.set_ylim((-1 , 30))
ax.set_xticks (np.linspace (0, 1, 6))
ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# ax.set_yticks (np.linspace (0.50,1,6))
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0f}'))
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
ax.yaxis.set_minor_locator(tck.AutoMinorLocator())

plt.savefig   ("m-x-contacts.png", dpi=1000)
