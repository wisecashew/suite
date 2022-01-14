'''
plot the energy dump file created by Monte Carlo Engine
''' 
 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="Read an energy dump file and calculate heat capacity.")
parser.add_argument('-i', metavar=': input energy dump file (energydump.txt)', dest='i', action='store', help='enter address of coordinate file')
parser.add_argument('-eh', metavar=': name of energy histogram.', type=str, dest='e', action='store', help='enter the edge length of the cubic simulation box')
parser.add_argument('-c', metavar=': text file containing literally just the heat capacity.', type=str, dest='c', action='store', help='file name of heat capacity')
parser.add_argument('-T', metavar=': thermodynamic temperature', type=float, dest='T', action='store', help='enter the thermodynamic temperature T of the simulation with k = 1 ')
args = parser.parse_args() 

with open("energydump.txt", 'r') as f:
    for line in f:
        if line.startswith("Energy"):
            break
    df = pd.read_csv(f , sep='|', names =["Energy", "Step_Number"] )
    f.close() 
    
energy = df["Energy"].values
step_num = df["Step_Number"].values 

plt.hist(energy, density=True)
plt.xlabel('Energy of system (in kT)')
plt.ylabel('Frequency')
plt.savefig(args.e, dpi=1200)
# plt.show( )

energy2 = energy**2 

Cv = ((np.mean(energy2) - (np.mean(energy))**2))/(args.T)**2 

f = open(args.c, 'w')
f.write("heat capacity - {}".format(Cv))
f.close()

print("heat capacity at constant volume is {:.2f}".format(Cv) )