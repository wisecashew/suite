#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import os 
import aux 


''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''


parser = argparse.ArgumentParser(description="Get the heat capacity for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting point for plotting.', default=100)

args = parser.parse_args()  

U_list = aux.dir2U ( os.listdir(".") ) 

dop = args.dop

# plt.figure(1)

i=1
max_contacts = (0,0,0)
for U in U_list:
    
    # print ("Currently plotting out stuff in U = " + str(U) + "...", end =' ', flush=True )

    temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(args.dop) ) ) 
    # print (temperatures)
    for temp in temperatures: 
        
        num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(args.dop) + "/" + str(temp) ) ) )
        # print(num_list) 
        for num in num_list:
            df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/"+"energydump"+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot","ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)
            # print (str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/"+args.e+"_"+str(num))
            if np.max(df["mm_tot"].values) > max_contacts[0]:
                max_contacts = ( np.max( df["mm_tot"].values ), U, temp) 

    
print ("N = " + str(dop) + ". The tuple is: ", max_contacts, flush=True)    
