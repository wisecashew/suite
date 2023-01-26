#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import os 
import aux 


''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''


parser = argparse.ArgumentParser(description="Get the heat capacity for simulation for every energy surface.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide a potential energy surface.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting point for plotting.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt') 
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.',default=False) 

args = parser.parse_args()  

dop_list = aux.dir2dop ( os.listdir(args.U) ) 

U = args.U 

dop_list = dop_list[0:7]
ldop_list = len(dop_list)

fig = plt.figure()
ax = plt.axes()

i=1
for dop in dop_list:
    
    print ("Currently plotting out stuff in N = " + str(dop) + "...", end =' ', flush=True )
    energy_list = np.asarray([])
    cv_list = np.asarray([])
    cv_err  = np.asarray([])
    cv_mean = np.asarray([])

    temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(dop) ) ) 

    for temp in temperatures: 
    
        num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(dop) + "/" + str(temp) ) ) )
        # print(num_list) 
        for num in num_list:
            df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot","ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)
            energy_list = df["energy"].values# [args.s:] 
            cv = (np.mean( (energy_list) ** 2 ) - np.mean ( energy_list )**2)/((temp**2)*dop) 
            cv_list = np.hstack ( (cv_list, cv ) ) 
        
        cv_mean = np.hstack ( ( cv_mean, np.mean ( cv_list ) ) ) 
        cv_err  = np.hstack ( ( cv_err , np.std ( cv_list )/np.sqrt(100) ) ) 
        
    plt.errorbar ( temperatures, np.asarray ( cv_mean ), yerr = np.asarray (cv_err), fmt='o',\
            markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, color=cm.winter(i/ldop_list) ) 

    i += 1
    print ("done!", flush=True) 


legend =  ["N = " + str(j) for j in dop_list ] 
plt.legend ( legend )
ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(10) )
ax.yaxis.get_major_locator().set_params(integer=True)
plt.xscale ('log')
plt.ylabel("$C_v/N$", fontsize=18)
plt.xlabel("Temperature (reduced)", fontsize=18)
plt.savefig ( "U_"+str(U)+"_CV.png", dpi=800)

if args.sp:
    plt.show()



