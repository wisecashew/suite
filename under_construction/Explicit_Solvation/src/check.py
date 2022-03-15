import pandas as pd 
import numpy as np 

df = pd.read_csv('energydump.txt', engine='python', names=['energy','mc','ac','nc','ts'], sep=' \| ') 

for i in range(len(df["energy"].values)):
    e = df["energy"].values[i]
    m = df["mc"].values[i] 
    a = df["ac"].values[i] 
    n = df["nc"].values[i]
    if (e != -m  -a  -n):
        print("Fuck up at " + str(df["ts"].values[i])) 
        break
