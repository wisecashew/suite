#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np
import aux


x = aux.get_f_of_s ( "U1", 0.01, 1, 32, "coords", 180000000)

print ("persistence length is =", x[0])
print ("r2 = ",x[1])


