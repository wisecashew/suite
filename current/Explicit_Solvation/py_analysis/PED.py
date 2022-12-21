#!/usr/licensed/anaconda3/2020.7/bin/python
'''
plot the energy dump file created by Monte Carlo Engine
''' 
 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Read an energy dump file and calculate heat capacity.")
parser.add_argument('-i', metavar=': input energy dump file (energydump.txt)', dest='i', action='store', help='enter address of coordinate file')
# parser.add_argument('-eh', metavar=': name of energy histogram.', type=str, dest='e', action='store', help='enter the edge length of the cubic simulation box')
parser.add_argument('-c', metavar=': text file containing literally just the heat capacity.', type=str, dest='c', action='store', help='file name of heat capacity')
parser.add_argument('-T', metavar=': thermodynamic temperature', type=float, dest='T', action='store', help='enter the thermodynamic temperature T of the simulation with k = 1 ')
parser.add_argument('-s', metavar=': value after which energy will be considered', type=int, dest='s', action='store', help='enter the starting point of trajectory. Deafult = 100.', default=100)
args = parser.parse_args() 

df = pd.read_csv(args.i, sep='\|', engine='python', names =["Energy","mm_contacts", "a_contacts", "n_contacts", "Step_Number"] )
    
energy = df["Energy"].values[args.s:]
step_num = df["Step_Number"].values[args.s:]

plt.hist(energy, density=True)
plt.xlabel('Energy of system (in kT)')
plt.ylabel('Frequency')
# plt.savefig(args.e, dpi=1200)
# plt.show( )

energy2 = energy**2 

print("temperature is " + str(args.T))
# print(energy)
Cv = ((np.mean(energy2) - (np.mean(energy))**2))/(args.T)**2 

f = open(args.c, 'w')
f.write("heat capacity - {}".format(Cv))
f.close()

print("heat capacity at constant volume is {:.2f}".format(Cv) )
