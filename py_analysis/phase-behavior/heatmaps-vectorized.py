#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import numba as nb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import numexpr as ne


# @nb.njit # (parallel=True)
def zmm(emmn, emma, g, T):
    return (1-g)*np.exp(-emmn/T, dtype=np.float128) + g * np.exp(-emma/T, dtype=np.float128)
    # return ne.evaluate("(1-g) * exp(-emmn/T) + g * exp(-emma/T)")

# @nb.njit # (parallel=True)
def zms(emsn, emsa, g, T):
    return (1-g)*np.exp(-emsn/T, dtype=np.float128) + g * np.exp(-emsa/T, dtype=np.float128)
    # return ne.evaluate("(1-g) * exp(-emsn/T) + g * exp(-emsa/T)")

# @nb.njit # (parallel=True)
def fmma(emmn, emma, g, T):
    z = zmm(emmn, emma, g, T)
    return g*np.exp(-emma/T, dtype=np.float128) / z
    # return ne.evaluate("g * exp(-emma/T) / z")

# @nb.njit # (parallel=True)
def fmsa(emsn, emsa, g, T):
    z = zms(emsn, emsa, g, T)
    return g*np.exp(-emsa/T, dtype=np.float128) / z
    # return ne.evaluate("g * exp(-emsa/T) / z")

# @nb.njit # (parallel=True)
def chi(emmn, emma, emsn, emsa, g, pv, T):
    fmsa_val = fmsa(emsn, emsa, g, T)
    fmma_val = fmma(emmn, emma, g, T)
    ems = pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn
    emm = 0.5* ( pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn ) 
    return 24/T * (ems-emm) # ((pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn) - 0.5*(pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn))


if __name__=="__main__":

    fig, axes = plt.subplots(nrows=9, ncols=9, figsize=(8,6), constrained_layout=True)
    
    start = time.time()

    # define density of time points and range of plots
    linrange = 1000
    plot_lim = 10
    
    # define energies to plot things over 
    E_mm_a, E_ms_a = np.meshgrid (np.linspace(-plot_lim, plot_lim, linrange), np.linspace (-plot_lim, plot_lim, linrange))

    # get temperatures
    T  = np.logspace (-2, 2, 50)
    T_broadcast = np.broadcast_to (T, (E_mm_a.shape[0], E_mm_a.shape[1], len(T)))
    
    # get the other data structures
    g = 0.25
    G = g * np.ones (E_mm_a.shape)

    pv = 1.0
    PV = pv * np.ones (E_mm_a.shape)

    E_mm_a_expanded = E_mm_a[:, :, np.newaxis] 
    E_ms_a_expanded = E_ms_a[:, :, np.newaxis] 
    
    
    G_expanded      = G [:, :, np.newaxis]
    PV_expanded     = PV[:, :, np.newaxis]

    E_mm_n_list = [-1, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0]
    E_ms_n_list = [-1, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0]
    rcount = -1
    for e_mm_n in E_mm_n_list:

        rcount +=  1
        print (f"rcount = {rcount}")
        ccount  = -1
        E_mm_n = e_mm_n * np.ones (E_mm_a.shape)
        E_mm_n_expanded = E_mm_n[:, :, np.newaxis]

        for e_ms_n in E_ms_n_list:
            print (f"ccount = {ccount}")
            ccount += 1

            ax = axes[rcount][ccount]
            ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', labelsize=3, pad=2)
            ax.tick_params(axis='x', labelsize=3, width=2)
            ax.tick_params(axis='y', labelsize=3, width=2)


            E_ms_n = e_ms_n * np.ones (E_ms_a.shape)
            E_ms_n_expanded = E_ms_n[:, :, np.newaxis]


            print ("Calculating chis...")
            Z_expanded = chi (E_mm_n_expanded, E_mm_a_expanded, E_ms_n_expanded, E_ms_a_expanded, G_expanded, PV_expanded, T_broadcast)
            print ("Calculated!")

            print ("Processing the calculated chis...")
            Z_max = np.max(Z_expanded, axis=-1)
            Z_min = np.min(Z_expanded, axis=-1)
            hold  = (Z_max * Z_min < 0) & (Z_max > 0)
            Z     = np.where (hold, Z_max - Z_min, 0)
            print ("Processed!")

            print ("begin plotting...")
            try:
                ax.pcolormesh (E_ms_a, E_mm_a, Z, cmap='Reds', norm=colors.SymLogNorm(linthresh=0.001, vmin=0, vmax=3000), shading="auto")   
            except ValueError:
                ax.pcolormesh (E_ms_a, E_mm_a, Z, cmap='Reds', vmin=0, vmax=1)   

            ax.set_xlim (-plot_lim, plot_lim)
            ax.set_ylim (-plot_lim, plot_lim)
            ax.set_yticks([-plot_lim,0,plot_lim])
            ax.set_xticks([-plot_lim,0,plot_lim])
            ax.set_xticklabels (ax.get_xticks(), weight='bold')
            ax.set_yticklabels (ax.get_yticks(), weight='bold')
            ax.minorticks_on()
            print ("plotted!")

    fig.savefig ("fast-plots-longer.png", dpi=1200, bbox_inches="tight")
    
    stop = time.time()
    print(f"Time required by this computation is {stop-start} seconds.")


