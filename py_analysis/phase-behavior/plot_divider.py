#!/Users/satyendhamankar/opt/anaconda3/envs/GBCG/bin/python

import numpy as np
from scipy.stats import linregress
import pandas as pd

df = pd.read_csv ("boundaries_tophalf.csv", names=["Emma", "Emmn"], sep=',', engine='python', skiprows=1)
print (df)

x = df["Emmn"].values
y = df["Emma"].values 


X = np.array ([x,y])
L = len(x) // 4

r2_list = []
for i in range (4):
    m, b, r, p, se = linregress (X[0,i*L:(i+1)*L], X[1,i*L:(i+1)*L])
    r2 = r**2
    r2_list.append (r2)

X = X[:, L:2*L]
L = len(X[0,:]) // 4
r2_list = []
for i in range (4):
    m, b, r, p, se = linregress (X[0,i*L:(i+1)*L], X[1,i*L:(i+1)*L])
    r2 = r**2
    r2_list.append (r2)


print (f"r2_list = {r2_list}")
