#!/Users/satyendhamankar/opt/anaconda3/envs/CG/bin/python

import numpy as np
from scipy.stats import linregress
import copy
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv ("boundaries_tophalf.csv", names=["Emma", "Emmn"], sep=',', engine='python', skiprows=1)
print (df)

x = df["Emmn"].values
y = df["Emma"].values 


X_ori = np.array ([x,y])
X = copy.copy(np.array ([x,y]))

range_of_indices = [0, len(X[0])-1]


for j in range (2):
    L = (range_of_indices[1]-range_of_indices[0]) // 4
    r2_list = []
    for i in range (4):
        m, b, r, p, se = linregress (X[0,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L], X[1,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L])
        r2 = r**2
        r2_list.append (r2)
    print (f"j = {j}, r2_list = {r2_list}")
    elbow_loc = np.argmin (r2_list)
    print (f"elbow_loc = {elbow_loc}")
    # excise that portion from X to go deeper
    range_of_indices = [range_of_indices[0]+elbow_loc * L, range_of_indices[0]+(elbow_loc+1)*L] # np.arange (elbow_loc*L,(elbow_loc+1)*L)
    

m1, b1, r1, p1, se1 = linregress (X[0, 0:range_of_indices[0]], X[1, 0:range_of_indices[0]])

print (f"Slope of top half = {m1}")
print (f"R squared = {r1}")

m2, b2, r2, p2, se2 = linregress (X[0, range_of_indices[1]:], X[1, range_of_indices[1]:])

print (f"Slope of bottom half = {m2}")
print (f"R squared = {r2}")


df = pd.read_csv ("boundaries_rightarm.csv", names=["Emma", "Emmn"], sep=',', engine='python', skiprows=1)

x = df["Emmn"].values
y = df["Emma"].values 

m3, b3, r3, p3, se3 = linregress (x, y)
print (f"Slope of bottom half = {m3}")
print (f"R squared = {r3}")

plt.plot (X[0], X[1], color='coral')
plt.plot (X[0,range_of_indices[0]:range_of_indices[1]], X[1,range_of_indices[0]:range_of_indices[1]], color='steelblue')
plt.show()
print (range_of_indices)
