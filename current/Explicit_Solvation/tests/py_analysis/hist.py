#!/usr/licensed/anaconda3/2020.7/bin/python
'''
plot the energy dump file created by Monte Carlo Engine
''' 
 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Read some sort of data file with second column being the relvant property and generating a histogram.")
parser.add_argument('-i', metavar=': input energy dump file (energydump.txt)', dest='i', action='store', help='enter address of coordinate file')
parser.add_argument('-s', metavar=': number of configurations to discard', dest='s', action='store', help='enter start of histogram', type=int, default=0)
parser.add_argument('-vl', metavar=': throw in a vertical line', dest='vl', action='store', type=float, help='enter position where you want a vertical red line', default=0) 
parser.add_argument('-xlo', metavar=': define lower x limit of histogram', dest='xlo', action='store', type=float, help='enter position of lower x limit of histogram', default=0.7)
parser.add_argument('-xhi', metavar=': define upper x limit of histogram', dest='xhi', action='store', type=float, help='enter position of higher x limit of histogram', default=6.0)
parser.add_argument('-ylo', metavar=': define lower y limit of histogram', dest='ylo', action='store', type=float, help='enter position of lower y limit of histogram', default=-0.2)
parser.add_argument('-yhi', metavar=': define upper y limit of histogram', dest='yhi', action='store', type=float, help='enter position of higher y limit of histogram', default=5.0)
args = parser.parse_args() 

df = pd.read_csv(args.i, sep='\s+', engine='python', names =["TS", "Rg"] ) 

plt.hist( df["Rg"].values[int(args.s):], density=True, bins=np.arange(args.xlo,args.xhi,0.05) ) 
plt.xlim([args.xlo,args.xhi])
plt.xticks(np.arange(args.xlo,args.xhi,0.5))
plt.ylim([args.ylo,args.yhi])
plt.axvline(x=args.vl, color='r')
plt.show() 

