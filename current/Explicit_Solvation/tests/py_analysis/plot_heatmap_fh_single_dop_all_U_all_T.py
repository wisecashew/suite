#!/usr/licensed/anaconda3/2020.7/bin/python

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
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--excl-vol', dest='ev', action='store_true', help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on screen.', default=False)

args = parser.parse_args()   

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

@plt.FuncFormatter
def fake_log(x, pos):
    'The two args are the value and tick position'
    return r'$10^{%d}$' % (x)

if __name__=="__main__":
    
    # get the entire list of potential energy surfaces 
    U_list  = aux.dir2U ( os.listdir (".") )
    # U_list  = U_list[0:4]
    mm_matrix = np.zeros((len(U_list), 14))
    #instantiate plt figure
    plt.figure( figsize=(10,6) )
    ax = plt.axes() 
    ax.tick_params ( axis='x', labelsize=16, pad=10 )
    ax.tick_params ( axis='y', labelsize=16, pad=10 )

    # instantiate some pertinent variables
    i=0
    Tmax = []

    for U in U_list:
        print ("Currently plotting out stuff in U = " + str(U) + "...", end=' ', flush=True )
        mm_list = np.asarray([])
        mm_err  = np.asarray([])
        mm_mean = np.asarray([])
        temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(args.dop) ) )
        Tmax.append ( np.max(temperatures) )
        for temp in temperatures: 
            num_list = np.unique ( aux.dir2nsim (os.listdir ( str(U) + "/DOP_"+str(args.dop) + "/" + str(temp) ) ) )
            
            for num in num_list:
                df = pd.read_csv(str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine='python', skiprows=args.s)
                mm_list = np.hstack ( ( mm_list, ( df["mm_tot"].values ) - (args.dop-1) ) )

            
            mm_err  = np.hstack ( ( mm_err , np.std  ( mm_list ) / np.sqrt(40) ) ) 
            mm_mean = np.hstack ( ( mm_mean, np.mean ( mm_list ) ) )

        mm_matrix[i,:] = mm_mean 
        i += 1 
        print ("done!", flush=True)
   
    print (mm_matrix)
    plt.imshow(mm_matrix, cmap='coolwarm', interpolation='nearest', origin='lower',norm=MidpointNormalize(midpoint=80), extent=[0.1,100,-0.2,0.2], aspect='auto')
    my_cmap = matplotlib.cm.get_cmap ('coolwarm_r')
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1))
    
    ax = plt.axes()

    plt.ylabel("$\chi$", fontsize=18)
    plt.xlabel("Temperature (reduced)", fontsize=18)
    # plt.xticks( temperatures, fontsize=16 )
    cbar = plt.colorbar(sm, orientation='vertical')
    cbar.set_ticks( [0, 0.5, 1] )
    cbar.set_ticklabels( ["Coil", "Random", "Globule"] )
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ( "State of polymer chain", fontsize=18, rotation=270, labelpad=25 ) 
    ax.set_yticks (np.linspace(-0.2,0.2,5))
    ax.set_xticklabels (temperatures)
    ax.set_xscale('log')
    plt.savefig("DOP_"+str(args.dop)+"_heatmap.png", dpi=1000)
    plt.tight_layout()
    if args.sp:
        plt.show()
