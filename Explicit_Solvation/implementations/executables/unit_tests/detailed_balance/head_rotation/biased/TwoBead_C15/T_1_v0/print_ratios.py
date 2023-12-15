import numpy as np 
import re
import os

temp_str = "Temperature"
hits = "hits" 
energy = "Energy of system"

hits_no = []
hits_on = []

l_dirs = os.listdir(".")

count_n_to_o = sum(1 for s in l_dirs if s.startswith("n_to_o"))

for k in range(1,count_n_to_o+1):
    f = open("n_to_o_"+str(k), 'r') 
    for line in f:
        if re.search (temp_str, line):
            res = re.search ("\d+\.\d+|\d+", line)
            T = float(res.group(0))
        elif re.search(energy, line):
            res = re.search("\d+\.\d+|\d+", line)
            En  = float(res.group(0)) 
        elif re.search (hits, line):
            res = re.search ("\d+", line) 
            hits_no.append(float(res.group(0)))
            break 
    f.close() 
    f = open("o_to_n_"+str(k), 'r')
    for line in f:
        if re.search (energy, line):
            res = re.search("\d+\.\d+|\d+", line)
            Eo  = float(res.group(0)) 
        elif re.search (hits, line):
            res = re.search ("\d+", line)
            hits_on.append(float(res.group(0)))
            break 

    f.close() 
truth = np.exp(1/T*(En-Eo))
pred = np.array(hits_on)/np.array(hits_no)
print (f"Ratio of empirical transitions (predicted) = {pred}")
print (f"Ratio of unbiased transitions (truth) = {np.exp(1/T*(En-Eo))}")
print (f"Mean transition rate = {np.mean(pred)}")
print (f"Standard error = {np.std(pred)/np.sqrt(len(pred))}")
print ("Error% = {}".format(np.abs(np.mean(pred)-truth)/truth*100))

