#!/home/satyend/.conda/envs/data_analysis/bin/python

import pandas as pd 
import numpy as np 
import matplotlib
import matplotlib.ticker as mticker
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
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)

args = parser.parse_args()   

divnorm = matplotlib.colors.TwoSlopeNorm(vmin=0, vcenter=60, vmax=350)

class MidpointNormalize(matplotlib.colors.Normalize):
    ## class from the mpl docs:
    # https://matplotlib.org/users/colormapnorms.html

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))    

def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"  # remove int() if you don't use MaxNLocator
    # return f"{10**val:.2e}"      # e-Notation



if __name__=="__main__":
    
    # get the entire list of potential energy surfaces 

    #instantiate plt figure
    fig = plt.figure( figsize=(8,6) )
    ax = fig.add_subplot(projection='3d' ) 
    ax.tick_params ( axis='x', labelsize=12 )
    ax.tick_params ( axis='y', labelsize=12 )
    ax.tick_params ( axis='z', labelsize=12 ) 

    # instantiate some pertinent variables
    i=0
    Tmax = []
    
    set_list = aux.dir2set (os.listdir (".") )
    # set_list = [set_list[2]]
    print (set_list)
    for s in set_list:
        mm_max   = -1
        chi_a_list = []
        chi_n_list = []
        print ("Processing stuff in set" + str(s) + "...", flush=True)
        U_list = aux.dir2U ( os.listdir ("set"+str(s)) )
        for U in U_list:
            if U == "U5" or U == "U7":
                continue
            print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True )
            mm_list = np.asarray([])
            mm_err  = np.asarray([])
            mm_mean = np.asarray([])
            temperatures = aux.dir2float ( os.listdir("set" + str(s) + "/" + str(U) +"/DOP_"+str(args.dop) ) )
            Tmax.append ( np.max(temperatures) )
            for temp in temperatures: 
                num_list = np.unique ( aux.dir2nsim (os.listdir ("set" + str(s) + "/" + str(U) + "/DOP_"+str(args.dop) + "/" + str(temp) ) ) )
                
                for num in num_list:
                    df = pd.read_csv("set" + str(s) + "/" + str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)
                    mm_list = np.hstack ( ( mm_list, ( df["mm_tot"].values ) - (args.dop-1) ) )

                # if U == "U1" and temp == 0.1:
                #     mm_max = np.max( mm_list )
                mm_err  = np.hstack ( ( mm_err , np.std  ( mm_list ) / np.sqrt(40) ) ) 
                mm_mean = np.hstack ( ( mm_mean, np.mean ( mm_list ) ) )
            
            chi_n = aux.get_chi_entropy("set" + str(s) + "/" + str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/geom_and_esurf.txt")[1] 
            chi_a = aux.get_chi_entropy("set" + str(s) + "/" + str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/geom_and_esurf.txt")[0] 
            chi_a_list = np.ones (len(temperatures))*chi_a
            chi_n_list = np.ones (len(temperatures))*chi_n
            ax.scatter3D(chi_a_list[::2], chi_n_list[::2], np.log( temperatures)[::2] , marker='o', edgecolors='k', c=mm_mean[::2], cmap='coolwarm_r', label='_nolegend_', norm=divnorm, depthshade=False)   
            print (mm_mean)
            i += 1
            # mm_max.append( np.max(mm_list) ) 
            print ("done!", flush=True)
            print ( "chi_a_list is ", chi_a_list )
            print ( "chi_n_list is ", chi_n_list )
    
    my_cmap = cm.coolwarm
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1))

    ax.set_ylabel("$\chi ^{a'}$", fontsize=14, labelpad=12)
    ax.set_xlabel("$\chi ^a$", fontsize=14, labelpad=12)
    ax.set_zlabel("$T^*$", fontsize=14, labelpad=12)
    
    cbar = plt.colorbar(sm, orientation='vertical', pad=0.2)
    cbar.set_ticks( [0, 0.5, 1] )
    cbar.set_ticklabels( ["Globule", "Random", "Coil"] )
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ( "State of polymer chain", fontsize=18, rotation=270, labelpad=25 ) 
    ax.grid(linewidth=20)
    ax.set_zticklabels ([0.01, 0.1, 0.5, 1, 5, 10, 50, 100])
    plt.savefig("DOP_"+str(args.dop)+"_scatter3d_mmcorr.png", dpi=1000)

    if args.sp:
        plt.show()

