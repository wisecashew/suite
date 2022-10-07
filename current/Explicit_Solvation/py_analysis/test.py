import numpy as np
import aux
import scipy.spatial.distance as ssd 

# glob = np.asarray([ [9,0,2], [8,0,1], [8,0,2], [8,1,1], [9,0,1], [9,1,2], [9,1,1], [8,1,2] ])
# coil = np.asarray([ [1,8,4], [0,9,5], [1,0,4], [2,1,3], [2,1,2], [3,1,1], [2,2,0], [1,1,9] ]) 
# coil = np.asarray([ [0,0,0], [0,0,1], [0,0,2], [0,0,3], [0,0,4], [0,0,5], [0,0,6], [0,0,7] ]) 

glob   = np.asarray ([[1,0,17],[1,17,0],[1,0,1], [0,1,0], [0,0,17], [0,0,0], [1,17,17], [2,0,0], [1,1,0], [1,0,0], [1,1,1],[2,1,0], [2,0,1], [2,17,0], [2,0,17],[1,1,17]])
coil   = np.asarray ( [[0,0,0],[0,0,1],[0,0,2],[0,0,3],[0,0,4],[0,0,5],[0,0,6],[0,0,7],[0,0,8],[0,0,9],[0,0,10],[0,0,11],[0,0,12],[0,0,13],[0,0,14],[0,0,15]] )

N = len(glob)

print ("length of glob is ", len(glob))
print ("length of coil is ", len(coil))

glob = aux.unfuck_polymer ( glob, 18, 18, 18 ) 
coil = aux.unfuck_polymer ( coil, 18, 18, 18 ) 

# print (glob)
# print (coil)

# calculate rg for globule 

r_com_glob = np.mean ( glob, axis=0 ) 
offset = glob - r_com_glob 
print ("glob com: ", r_com_glob)
rg_glob = np.sqrt ( np.sum (np.square (offset) / N ) ) 


# calculate rh for globule 
rh_glob = 1 / ( np.sum ( 1/ssd.pdist( glob, 'euclidean' ) )/ (N*(N-1)/2) )

print ("For the globular case...")
print ("Radius of gyration is {:.2f} and hydrodynamic radius is {:.2f}.".format (rg_glob, rh_glob) ) 
print ("Rg/Rh is ", rg_glob/rh_glob )


# calculate rg for globule 
r_com_coil = np.mean ( coil, axis=0 ) 
offset = coil - r_com_coil 
print ("coil com: ", r_com_coil)
rg_coil = np.sqrt ( np.sum (np.square (offset) / N ) ) 


# calculate rh for globule 
rh_coil = 1/ ( np.sum ( 1/ssd.pdist( coil, 'euclidean' ) )/ (N*(N-1)/2) ) 

print ("For the coil case...")
print ("Radius of gyration is {:.2f} and hydrodynamic radius is {:.2f}.".format (rg_coil, rh_coil) ) 
print ("Rg/Rh is ", rg_coil/rh_coil )


