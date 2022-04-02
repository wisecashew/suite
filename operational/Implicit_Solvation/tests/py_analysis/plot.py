#!/usr/licensed/anaconda3/2020.7/bin/python
'''
plot the energy dump file created by Monte Carlo Engine
''' 
 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Read a two column file and plotting the second column (y-axis) against the first column (x-axis).")
parser.add_argument('-i', metavar=': input energy dump file (energydump.txt)', dest='i', action='store', help='enter address of coordinate file')
parser.add_argument('-ylo', metavar=': lower ylimit on plot.', dest='ylo', type=float, action='store', help='enter lower limit of plot', default=-0.2)
parser.add_argument('-yhi', metavar=': higer ylimit on plot.', dest='yhi', type=float, action='store', help='enter lower limit of plot', default=7)
parser.add_argument('-xlo', metavar=': lower xlimit on plot.', dest='xlo', type=float, action='store', help='enter lower limit of plot', default=0)
parser.add_argument('-xhi', metavar=': higer xlimit on plot.', dest='xhi', type=float, action='store', help='enter lower limit of plot', default=7)
args = parser.parse_args() 

df = pd.read_csv(args.i, sep='\s+', engine='python', names =["Temp", "Cv"] ) 

plt.plot( df["Temp"].values, df["Cv"].values, linestyle='dashed', marker='^' ) 
plt.xlim([args.xlo, args.xhi])
plt.ylim([args.ylo, args.yhi])
plt.axhline(y=0, color='r')
plt.show() 

