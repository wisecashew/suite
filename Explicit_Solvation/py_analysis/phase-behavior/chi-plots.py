#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 



if __name__=="__main__":

    fig = plt.figure()
    ax  = plt.axes ()

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
    E_mm_n = np.arange (elow, ehigh+0.1, 0.1)
    T_range = np.logspace (-2, 2, 25)

    for e_mm_n in E_mm_n:
        rgba_color = cm.GnBu (norm(e_mm_n))
        plt.plot (T_range, chi(E_mm_a, e_mm_n, E_ms_a, E_ms_n, T_range), marker='o', markeredgecolor='k', c=rgba_color)
        # print (chi(E_mm_a, e_mm_n, E_ms_a, E_ms_n, T_range))
    ax.set_xscale ("log")
    ax.set_yscale ("symlog")
    ax.axhline (y=0, c='crimson')
    ax.minorticks_on()
    # ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.show()
