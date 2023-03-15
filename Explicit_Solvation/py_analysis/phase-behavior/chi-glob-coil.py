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

    fig = plt.figure(figsize=(6,4))
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=16)
    ax.tick_params(axis='x', labelsize=16, pad=10)
    ax.tick_params(axis='y', labelsize=16)

    ehigh  = -2
    elow   = -3
    g = 0.2514
    zmm  = lambda emma, emmn, T: g*np.exp (-1/T * emma) + (1-g)*np.exp (-1/T * emmn)
    zms  = lambda emsa, emsn, T: g*np.exp (-1/T * emsa) + (1-g)*np.exp (-1/T * emsn)
    fmma = lambda emma, emmn, T: g*np.exp (-1/T * emma)/zmm(emma, emmn, T)
    fmsa = lambda emsa, emsn, T: g*np.exp (-1/T * emsa)/zms(emsa, emsn, T)
    norm = matplotlib.colors.Normalize (vmin=elow, vmax=ehigh)
    def chi (emma, emmn, emsa, emsn, T):
        t1 = fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn
        t2 = fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn
        return 24/T * (t1 - 0.5 * t2)


    E_mm_a = -3; E_ms_a = -1.4; E_ms_n = -1.4;
    E_mm_n = np.arange (elow, ehigh+0.1, 0.2)
    T_range = np.logspace (-2, 2, 25)

    for e_mm_n in E_mm_n:
        rgba_color = cm.summer (norm(e_mm_n))
        plt.plot (T_range, chi(E_mm_a, e_mm_n, E_ms_a, E_ms_n, T_range), marker='o', markeredgecolor='k', c=rgba_color)
        

    ax.set_xscale ("log")
    ax.set_yscale ("symlog")
    ax.set_ylim (bottom=-10, top=500)
    ax.set_xlim (left=0.01,right=100)
    ax.axhline (y=0, c='crimson')
    fig.tight_layout()
    ax.set_xticks ([1e-2, 1e-1, 1, 1e+1, 1e+2])
    ax.set_yticklabels (["$\\mathbf{-10^{1} }$", "$\\mathbf{-10^{0} }$", "$\\mathbf{0}$", "$\\mathbf{10^{0} }$", "$\\mathbf{10^{1} }$", "$\\mathbf{10^{2} }$"], weight='bold')
    ax.set_xticklabels (["$\\mathbf{10^{-2} }$", "$\\mathbf{10^{-1} }$", "$\\mathbf{10^{0} }$", "$\\mathbf{10^{1} }$", "$\\mathbf{10^{2} }$"], weight='bold')
    ax.minorticks_on()

    yaxis = plt.gca().yaxis
    yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.0e}'))
    # plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
   
    plt.savefig ("g-to-c.png", bbox_inches='tight', dpi=1200)


