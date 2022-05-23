import aux 
import numpy as np

master_dict = {} 

v = np.asarray( [ [0, 1, 2], [0, 2, 2], [0, 2, 1], [0, 1, 1], [9, 1, 1], [9, 1, 2], [9, 2, 2], [9, 2, 1] ]  )
v = aux.unfuck_polymer (v, 10, 10,10)
master_dict[0] = [ v ]

hydr = aux.get_Rh ( master_dict, 10, 10, 10)

print ( "hydrodynamic radius is {:.2f}".format(hydr))
