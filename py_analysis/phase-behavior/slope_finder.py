#!/Users/satyendhamankar/opt/anaconda3/bin/python

import numpy as np
import pandas as pd
import copy 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import scipy.optimize as so
import time

if __name__=="__main__":

    df = pd.read_csv ("boundaries_top_left_half_c_to_g.csv", names=["Emma", "Emmn"], sep=',', engine='python', skiprows=1)

    x = df["Emmn"].values
    y = df["Emma"].values

    X = np.array ([x, y])
    range_of_indices = [0, len(X[0])-1]
    r2_list = []

    for j in range(2):
        L = (range_of_indices[1] - range_of_indices[0]) // 4
        r2_list.clear() 
        for i in range(4):
            m, b, r, p, se = linregress (X[0,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L], X[1,range_of_indices[0]+i*L:range_of_indices[0]+(i+1)*L])
            r2 = r**2
            r2_list.append (r2)

        elbow_loc = np.argmin (r2_list) 
        range_of_indices = [range_of_indices[0]+elbow_loc * L, range_of_indices[0]+(elbow_loc+1)*L] # np.arange (elbow_loc*L,(elbow_loc+1)*L)

    m1, b1, r1, p1, se1 = linregress (X[0, 0:range_of_indices[0]], X[1, 0:range_of_indices[0]])
    m2, b2, r2, p2, se2 = linregress (X[0, range_of_indices[1]:] , X[1, range_of_indices[1]:] )

    print (f"Slope of bottom    = {m1}")
    print (f"r squared          = {r1}")
    print (f"Slope of top       = {m2}")
    print (f"r squared          = {r2}")