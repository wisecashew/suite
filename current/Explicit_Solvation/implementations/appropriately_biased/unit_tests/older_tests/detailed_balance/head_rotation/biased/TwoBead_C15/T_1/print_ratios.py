import numpy as np 
import re

temp_str = "Temperature"
hits = "hits" 
energy = "Energy of system"

for k in ["1", "2", "3"]:
    f = open("n_to_o_"+k, 'r') 
    for line in f:
        if re.search (temp_str, line):
            res = re.search ("\d+\.\d+|\d+", line)
            T = float(res.group(0))
        elif re.search(energy, line):
            res = re.search("\d+\.\d+|\d+", line)
            En  = float(res.group(0)) 
        elif re.search (hits, line):
            res = re.search ("\d+", line) 
            hits_no = float(res.group(0)) 
            break 
    f.close() 
    f = open("o_to_n_"+k, 'r')
    for line in f:
        if re.search (energy, line):
            res = re.search("\d+\.\d+|\d+", line)
            Eo  = float(res.group(0)) 
        elif re.search (hits, line):
            res = re.search ("\d+", line)
            hits_on = float(res.group(0))
            break 

    f.close() 
    pred = hits_on/hits_no
    truth = np.exp(1/T*(En-Eo))
    print (k+". Ratio of empirical transitions (predicted) = ", hits_on/hits_no)
    print (k+". Ratio of unbiased transitions (truth) = ", np.exp(1/T*(En-Eo)))
    print ("Error% = {}".format(np.abs(pred-truth)/truth*100))

