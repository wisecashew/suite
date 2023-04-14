#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import multiprocessing as mp
import numexpr as ne

def zmm(emmn, emma, g, T):
    return ne.evaluate("(1-g) * exp(-emmn/T) + g * exp(-emma/T)")

def zms(emsn, emsa, g, T):
    return ne.evaluate("(1-g) * exp(-emsn/T) + g * exp(-emsa/T)")

def fmma(emmn, emma, g, T):
    z = zmm(emmn, emma, g, T)
    return ne.evaluate("g * exp(-emma/T) / z")

def fmsa(emsn, emsa, g, T):
    z = zms(emsn, emsa, g, T)
    return ne.evaluate("g * exp(-emsa/T) / z")

def chi(emmn, emma, emsn, emsa, g, pv, T):
    fmsa_val = fmsa(emsn, emsa, g, T)
    fmma_val = fmma(emmn, emma, g, T)
    ems = pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn
    emm = 0.5* ( pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn ) 
    return 24/T * (ems-emm) # ((pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn) - 0.5*(pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn))

def chi_min (emmn, emma, emsn, emsa, g, pv, T):
    fmsa_val = fmsa(emsn, emsa, g, T)
    fmma_val = fmma(emmn, emma, g, T)
    ems = pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn
    emm = 0.5* ( pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn ) 
    chi_val = 24/T * (ems-emm)
    chi_min = np.min (chi_val)
    return chi_min

def chi_max (emmn, emma, emsn, emsa, g, pv, T):
    fmsa_val = fmsa(emsn, emsa, g, T)
    fmma_val = fmma(emmn, emma, g, T)
    ems = pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn
    emm = 0.5* ( pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn ) 
    chi_val = 24/T * (ems-emm)
    chi_max = np.max (chi_val)
    return chi_max




def make_heatmap_v1(E_mm_n, E_ms_n, T_range, linrange, plot_lim, g, pv):

    y, x = np.meshgrid(np.linspace(-plot_lim, plot_lim, linrange), np.linspace(-plot_lim, plot_lim, linrange))
    z = np.zeros(x.shape)
    for x_idx in range(linrange):
        for y_idx in range(linrange):
            chi_T = chi(E_mm_n, y[x_idx, y_idx], E_ms_n, x[x_idx, y_idx], g, pv, T_range)
            chi_min = np.min(chi_T)
            chi_max = np.max(chi_T)
            if (chi_min * chi_max) < 0 and chi_max > 0:
                z[x_idx, y_idx] = chi_max - chi_min
            elif (chi_min * chi_max) > 0:
                z[x_idx, y_idx] = 0
            elif chi_min * chi_max == 0:
                z[x_idx, y_idx] = 0
            else:
                z[x_idx, y_idx] = 0
    z_min = np.min(z)
    z_max = np.max(z)
    

    fig = plt.figure(figsize=(4/1.6, 3/1.6), constrained_layout=True)
    ax  = plt.axes()

    try:
        ax.pcolormesh (x, y, z, cmap='Reds', norm=colors.SymLogNorm(linthresh=0.001, vmin=0, vmax=3000), shading="auto")   
    except ValueError:
        ax.pcolormesh (x, y, z, cmap='Reds', vmin=0, vmax=1)  


    fig.savefig ("fast-plots", dpi=1200, bbox_inches="tight")
    return


def make_heatmap_v2(E_mm_n, E_ms_n, T_range, linrange, plot_lim, g, pv):
    
    y, x = np.meshgrid(np.linspace(-plot_lim, plot_lim, linrange), np.linspace(-plot_lim, plot_lim, linrange))

    T_tile = np.tile(T_range, (linrange, linrange, 1)).transpose(2, 0, 1)
    emma_tile = np.tile(y, (linrange, linrange, 1))
    emmn_tile = np.tile(E_mm_n, (linrange, linrange, len(T_range)))
    emsa_tile = np.tile(x, (linrange, linrange, 1))
    emsn_tile = np.tile(E_ms_n, (linrange, linrange, len(T_range)))
    g_tile = np.tile(g, (linrange, linrange, len(T_range)))
    pv_tile = np.tile(pv, (linrange, linrange, len(T_range)))

    z = chi(emmn_tile, emma_tile, emsn_tile, emsa_tile, g_tile, pv_tile, T_tile)
    
    chi_min = np.min(z, axis=2)
    chi_max = np.max(z, axis=2)
    z[(chi_min * chi_max) < 0] = chi_max[(chi_min * chi_max) < 0] - chi_min[(chi_min * chi_max) < 0]
    z[(chi_min * chi_max) >= 0] = 0
    
    z_min = np.min(z)
    z_max = np.max(z)

    fig = plt.figure(figsize=(4/1.6, 3/1.6), constrained_layout=True)
    ax  = plt.axes()

    try:
        ax.pcolormesh (x, y, z, cmap='Reds', norm=colors.SymLogNorm(linthresh=0.001, vmin=0, vmax=3000), shading="auto")   
    except ValueError:
        ax.pcolormesh (x, y, z, cmap='Reds', vmin=0, vmax=1)  

    fig.savefig ("fast-plots", dpi=1200, bbox_inches="tight")
    return

def make_heatmap_v3(E_mm_n, E_ms_n, T_range, linrange, plot_lim, g, pv):
    
    y, x = np.meshgrid(np.linspace(-plot_lim, plot_lim, linrange), np.linspace(-plot_lim, plot_lim, linrange))

    # y = E_mm_a
    # x = E_ms_a 

    E_mm_n = E_mm_n * np.ones (x.shape)
    E_ms_n = E_ms_n * np.ones (x.shape)
    g      = g * np.ones  (x.shape)
    pv     = pv * np.ones (x.shape)

    z_min = chi_min(E_mm_n, y, E_ms_n, x, g, pb, T_range) # emsn_tile, emsa_tile, g_tile, pv_tile, T_tile)
    
    chi_min = np.min(z, axis=2)
    chi_max = np.max(z, axis=2)
    z[(chi_min * chi_max) < 0] = chi_max[(chi_min * chi_max) < 0] - chi_min[(chi_min * chi_max) < 0]
    z[(chi_min * chi_max) >= 0] = 0
    
    z_min = np.min(z)
    z_max = np.max(z)

    fig = plt.figure(figsize=(4/1.6, 3/1.6), constrained_layout=True)
    ax  = plt.axes()

    try:
        ax.pcolormesh (x, y, z, cmap='Reds', norm=colors.SymLogNorm(linthresh=0.001, vmin=0, vmax=3000), shading="auto")   
    except ValueError:
        ax.pcolormesh (x, y, z, cmap='Reds', vmin=0, vmax=1)  

    fig.savefig ("fast-plots", dpi=1200, bbox_inches="tight")
    return



if __name__=="__main__":

    start = time.time()
    g  = 0.25
    pv = 1.0
    linrange = 100
    T_range  = np.logspace (-2, 2, linrange)
    plot_lim = 1
    E_mm_n = [-1,0,1]
    E_ms_n = [-1,0,1]

    y, x = np.meshgrid (np.linspace (-plot_lim,plot_lim,linrange), np.linspace (-plot_lim,plot_lim,linrange) )

    make_heatmap_v2 (E_mm_n, E_ms_n, T_range, linrange, plot_lim, g, pv)
    stop = time.time()

    print(f"Timerequired by this computation is {stop-start} seconds.")


