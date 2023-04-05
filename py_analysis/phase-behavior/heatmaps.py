#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator
import matplotlib.colors as colors

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
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

 
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
        return 24 / T* (t1-0.5 * t2) # 24/T * (t1 - 0.5 * t2)  

    def make_heatmap(x, y, T, g, pv, linrange):

        
        E_mm_a = -3 
        E_mm_n = -3

        z = np.zeros (x.shape)

        for x_idx in range (linrange):
            for y_idx in range(linrange):
                chi_T = chi (E_mm_a, E_mm_n, x[x_idx, y_idx], y[x_idx, y_idx], g, pv, T_range)
                chi_min = np.min(chi_T) # * T_range[np.argmin(chi_T)]
                chi_max = np.max(chi_T) # * T_range[np.argmax(chi_T)]
                z[x_idx, y_idx] = chi_min * chi_max 

        return z

    g = 0.0
    pv = 1.0
    T_range  = np.logspace (-2,np.log10(25),100)
    linrange = 100
    plot_lim = 50
    y, x = np.meshgrid (np.linspace (-plot_lim,plot_lim*2,linrange), np.linspace (-plot_lim,plot_lim,linrange))


    z = make_heatmap(x, y, T_range, g, pv, linrange)

    z[z>0] = z[z>0]/np.max(z[z>0])
    z[z>0.8] = 1
    z_max = np.max(z)

    try:
        z[z<0] = -z[z<0]/np.min(z[z<0])
        z[z<-0.8] = -1
        z_min = np.min(z)
    except ValueError:
        z_min = -1
        pass
      

    print ("z_min = {}".format(z_min))
    print ("z_max = {}".format(z_max))
    ax.pcolormesh (x, y, z, cmap='RdBu', norm=colors.TwoSlopeNorm (vmin=z_min, vcenter=0, vmax=z_max) )  
    

    # ax.set_xscale ("log")
    # ax.set_yscale ("symlog")
    # ax.set_ylim (bottom=-500, top=500)
    # ax.set_xlim (left=0.01,right=100)
    # ax.axhline (y=0, c='crimson')
    # ax.set_xticks ([1e-2, 1e-1, 1, 1e+1, 1e+2])

    ax.minorticks_on()
    # ax.legend(fontsize="xx-small", loc="lower right")
   
    plt.savefig (f"chi-heatmap-pv-{pv}.png", bbox_inches='tight', dpi=1200)

        


