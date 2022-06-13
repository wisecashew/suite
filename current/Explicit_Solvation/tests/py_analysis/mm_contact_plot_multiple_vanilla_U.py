#!/usr/bin/env python3

import pandas as pd 
import numpy as np 
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import argparse 
import aux
import os 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

parser = argparse.ArgumentParser(description="Get the contacts for simulation for every energy surface, provided you give the volume fraction.")
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide a potential energy surface.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
# parser.add_argument('--excl-vol', dest='ev', action='store_true', help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)

args = parser.parse_args()   

if __name__=="__main__":
    
    # get the entire list of potential energy surfaces 
    U = args.U

    #instantiate plt figure
    plt.figure( figsize=(8,6) )
    
    # if we want the Uexcl contacts
    dop_list = aux.dir2dop ( os.listdir( U ) )

    # instantiate some pertinent variables
    i=0
    Tmax = []
    mm_max = [] 
    for dop in dop_list[0:6]:
        print ("Currently plotting out stuff in N = " + str(dop) + "...", end='\n', flush=True )
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])

        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
            num_list = np.unique ( aux.dir2nsim (os.listdir ( str(U) + "/DOP_"+str(dop) + "/" + str(temp) ) ) )
            
            for num in num_list:
                df = pd.read_csv(str(U)+"/DOP_"+str(dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)
                mm_list = np.hstack ( ( mm_list, ( df["mm_tot"].values ) - (dop-1) ) )
        
            mm_err  = np.hstack ( ( mm_err , np.std  ( mm_list ) / np.sqrt(5) ) ) 
            mm_mean = np.hstack ( ( mm_mean, np.mean ( mm_list ) ) )

        plt.errorbar(temperatures, np.asarray(mm_mean), yerr=np.asarray(mm_err), fmt='o', markeredgecolor='k', linestyle='-', elinewidth=1, capsize=0, color=cm.copper(i/9), label="N="+str(dop) )   
        i += 1
        mm_max.append( np.max(mm_list) ) 

        print ("done!", flush=True)
    
    ax = plt.axes() 
    plt.ylabel("$\\frac \\langle M_C \\rangle$", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    plt.xticks(np.arange(0.01, np.max(Tmax)+2))
    ytick_list = np.unique ( np.linspace (0, np.max(mm_max)+1, 20, dtype=int) )

    ytick_list = np.arange ( 0, ytick_list[-1], np.ceil(int(ytick_list[-1])/10 ) )  

    if not ( (np.max(mm_max)+1) in ytick_list):
        ytick_list = np.hstack ( (ytick_list, np.max(mm_max)+1 ) ) 

    # plt.yticks( ytick_list ) 
    plt.legend( [ str(d) for d in dop_list] )
    ax.set_xscale('log')
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12)
    plt.savefig(str(args.U)+"_mmcorr.png", dpi=800)

    if args.sp:
        plt.show()

