#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator

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

    g  = 0.25
    
    zmm  = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float128) + (1-g)*np.exp ((-1/T * emmn), dtype = np.float128)
    zms  = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float128) + (1-g)*np.exp ((-1/T * emsn), dtype = np.float128)
    fmma = lambda emma, emmn, T: g*np.exp ((-1/T * emma), dtype=np.float128) / zmm(emma, emmn, T)
    fmsa = lambda emsa, emsn, T: g*np.exp ((-1/T * emsa), dtype=np.float128) / zms(emsa, emsn, T)

    def chi (emma, emmn, emsa, emsn, pv, T):
        t1 = pv*(fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn) + (1-pv)*emmn
        return 24/T * (t1 - 0.5 * t2)   

    norm = matplotlib.colors.TwoSlopeNorm (vmin=0, vcenter=0.5, vmax=1.0)

    ems_list = [0]
    E_mm_a = 1
    E_mm_n = 0
    E_ms_n = 0

    T_range = np.logspace (-2, 2, 100)

    it = 0
    for pv in [1.0]:

        print (f"In pv = {pv}...")
        fig = plt.figure(it) #figsize=(4/1.6,3/1.6), constrained_layout=True)
        ax  = plt.axes ()
        ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=8)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        it += 1

        for E_ms_a in ems_list:
            rgba_color = cm.PiYG (norm(E_ms_a))
            c = chi(E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, T_range)
            print (f"c[0] = {c[0]}")
            if c[0] > 0:
                print ("Collapsed at cold temperature.")
                print ("max chi = {}".format(np.max(c)-c[0]))
                ax.plot (T_range, chi(E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, T_range) - chi(E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, 0.01), marker='o', markeredgecolor='k', c=rgba_color, label=f"$p_v = {pv}, E_{{ms}}^{{a}} = {E_ms_a}$")
            else:
                ax.plot (T_range, chi(E_mm_a, E_mm_n, E_ms_a, E_ms_n, pv, T_range), marker='o', markeredgecolor='k', c=rgba_color, label=f"$p_v = {pv}, E_{{ms}}^{{a}} = {E_ms_a}$")

        ax.set_xscale ("log")
        ax.set_yscale ("symlog")
        # ax.set_ylim (bottom=-500, top=500)
        ax.set_xlim (left=0.01,right=100)
        ax.axhline (y=0, c='crimson')
        ax.set_xticks ([1e-2, 1e-1, 1, 1e+1, 1e+2])
        ax.minorticks_on()
        ax.legend(fontsize="xx-small", loc="lower right")
        plt.savefig (f"dw-chi_pv_{pv}.png", bbox_inches='tight', dpi=1200)

        


