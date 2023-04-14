#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator
import matplotlib.colors as colors
import time

import argparse 
parser = argparse.ArgumentParser (description="Create images based on emsa and emsn.")
# parser.add_argument ("-g", dest='g', action='store', type=float, help="Enter g value.")
# parser.add_argument ("-p", dest='p', action='store', type=float, help="Enter p value.")

args = parser.parse_args()


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))



if __name__=="__main__":
    

    fig = plt.figure(figsize=(4/1.6,3/1.6), constrained_layout=True)
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=11)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

     
    zmm  = lambda emma, emmn, g, T: g*np.exp ((-1/T * emma), dtype=np.float128) + (1-g)*np.exp ((-1/T * emmn), dtype = np.float128)
    zms  = lambda emsa, emsn, g, T: g*np.exp ((-1/T * emsa), dtype=np.float128) + (1-g)*np.exp ((-1/T * emsn), dtype = np.float128)
    fmma = lambda emma, emmn, g, T: g*np.exp ((-1/T * emma), dtype=np.float128) / zmm(emma, emmn, g, T)
    fmsa = lambda emsa, emsn, g, T: g*np.exp ((-1/T * emsa), dtype=np.float128) / zms(emsa, emsn, g, T)
    chi  = lambda emma, emmn, emsa, emsn, g, pv, T: 24/T * ( (pv*(fmsa(emsa, emsn, g, T)*emsa + (1-fmsa(emsa, emsn, g, T))*emsn) + (1-pv)*emsn) - 0.5 *(pv*(fmma(emma, emmn, g, T)*emma + (1-fmma(emma, emmn, g, T))*emmn) + (1-pv)*emmn) )


    # def make_heatmap(x, y, T, g, pv, linrange):
    start = time.time()
    
    E_mm_n = -3 
    E_mm_a = -3
    E_ms_a = -25
    E_ms_n = 0
    # g  = 0.25
    # pv = 1.0
    T_range  = np.logspace (-2, 2, 100)
    linrange = 100
    plot_lim = 10

    y, x = np.meshgrid (np.linspace (0,1,linrange), np.linspace (0,1,linrange) )
    z = np.zeros (x.shape)

    for x_idx in range (linrange):
        for y_idx in range(linrange):
            chi_T = chi (E_mm_a, E_mm_n, E_ms_a, E_ms_n, x[x_idx,y_idx], y[x_idx, y_idx], T_range)
            chi_min = np.min(chi_T) # * T_range[np.argmin(chi_T)]
            chi_max = np.max(chi_T) # * T_range[np.argmax(chi_T)]
            if (chi_min * chi_max) < 0 and chi_max > 0:
                z[x_idx, y_idx] = np.max (chi_T)
                print ("pv = ",y[x_idx, y_idx])
                print ("g = ",x[x_idx, y_idx])
                print ("max chi_T  = ",np.max(chi_T))
                print ("min chi_T  = ",np.min(chi_T))
            elif (chi_min * chi_max) > 0:
                z[x_idx, y_idx] = 0
            elif chi_min * chi_max == 0:
                z[x_idx, y_idx] = 0
            else:
                z[x_idx, y_idx] = 0

    z_min = np.min(z)
    z_max = np.max(z)

    print ("z_min = {}".format(z_min))
    print ("z_max = {}".format(z_max))
    try:
        ax.pcolormesh (x, y, z, cmap='Reds', norm=colors.Normalize(vmin=z_min, vmax=z_max), shading="auto")   
    except ValueError:
        ax.pcolormesh (x, y, z, cmap='Reds', vmin=0, vmax=1)   

    ax.set_xlim (0, 1)
    ax.set_ylim (0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels (ax.get_xticks(), weight='bold')
    ax.set_yticklabels (ax.get_yticks(), weight='bold')
    ax.minorticks_on()
    
   
    plt.savefig (f"phases-g-pv.png", bbox_inches='tight', dpi=1200)
    stop = time.time()

    print ("Time for simulation {} seconds.".format (stop-start))
        
    """
    def zmm (emma, emmn, g, T):
        return g*np.exp ((-1/T * emma), dtype=np.float128) + (1-g)*np.exp ((-1/T * emmn), dtype = np.float128)

    def zms (emsa, emsn, g, T):
        return g*np.exp ((-1/T * emsa), dtype=np.float128) + (1-g)*np.exp ((-1/T * emsn), dtype = np.float128)

    def fmma (emma, emmn, g, T):
        return g*np.exp ((-1/T * emma), dtype=np.float128) / zmm(emma, emmn, g, T)

    def fmsa (emsa, emsn, g, T):
        return g*np.exp ((-1/T * emsa), dtype=np.float128) / zms(emsa, emsn, g, T)


    def chi (emma, emmn, emsa, emsn, g, pv, T):
        t1 = pv*(fmsa(emsa, emsn, g, T)*emsa + (1-fmsa(emsa, emsn, g, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, g, T)*emma + (1-fmma(emma, emmn, g, T))*emmn) + (1-pv)*emmn
        return 24 / T * (t1 - 0.5 * t2) # 24/T * (t1 - 0.5 * t2)  
    """

