#!~/.conda/envs/data_analysis/bin/python

import numpy as np
import numba as nb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import numexpr as ne


@nb.njit # (parallel=True)
def zmm(emmn, emma, g, T):
    return (1-g)*np.exp(-emmn/T) + g * np.exp(-emma/T)
    # return ne.evaluate("(1-g) * exp(-emmn/T) + g * exp(-emma/T)")

@nb.njit # (parallel=True)
def zms(emsn, emsa, g, T):
    return (1-g)*np.exp(-emsn/T) + g * np.exp(-emsa/T)
    # return ne.evaluate("(1-g) * exp(-emsn/T) + g * exp(-emsa/T)")

@nb.njit # (parallel=True)
def fmma(emmn, emma, g, T):
    z = zmm(emmn, emma, g, T)
    return g*np.exp(-emma/T) / z
    # return ne.evaluate("g * exp(-emma/T) / z")

@nb.njit # (parallel=True)
def fmsa(emsn, emsa, g, T):
    z = zms(emsn, emsa, g, T)
    return g*np.exp(-emsa/T) / z
    # return ne.evaluate("g * exp(-emsa/T) / z")

@nb.njit # (parallel=True)
def chi(emmn, emma, emsn, emsa, g, pv, T):
    fmsa_val = fmsa(emsn, emsa, g, T)
    fmma_val = fmma(emmn, emma, g, T)
    ems = pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn
    emm = 0.5* ( pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn ) 
    return 24/T * (ems-emm) # ((pv*(fmsa_val*emsa + (1-fmsa_val)*emsn) + (1-pv)*emsn) - 0.5*(pv*(fmma_val*emma + (1-fmma_val)*emmn) + (1-pv)*emmn))


if __name__=="__main__":

    E_ms_a_list = [-1, -0.5, 0] # [-1, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0]
    E_ms_n_list = [-1, -0.5, 0] # [-1, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0]

    fig, axes = plt.subplots(nrows=len(E_ms_a_list), ncols=len(E_ms_n_list), figsize=(8,6), constrained_layout=True)
    start = time.time()

    # define density of time points and range of plots
    linrange = 1000
    plot_lim = 5/2

    # define energies to plot things over 
    E_mm_a, E_mm_n = np.meshgrid (np.linspace(-plot_lim, plot_lim, linrange), np.linspace (-plot_lim, plot_lim, linrange))

    # get temperatures
    T  = np.logspace (-2, 2, 50)
    T_broadcast = np.broadcast_to (T, (E_mm_a.shape[0], E_mm_a.shape[1], len(T)))

    # get the other data structures
    g = 0.25
    G = g * np.ones (E_mm_a.shape)

    pv = 1.0
    PV = pv * np.ones (E_mm_a.shape)

    E_mm_a_expanded = E_mm_a[:, :, np.newaxis] 
    E_mm_n_expanded = E_mm_n[:, :, np.newaxis] 

    G_expanded      = G [:, :, np.newaxis]
    PV_expanded     = PV[:, :, np.newaxis]


    rcount = -1
    for e_ms_a in E_ms_a_list:

        rcount +=  1
        print (f"rcount = {rcount}", flush=True)
        ccount  = -1
        E_ms_a = e_ms_a * np.ones (E_mm_a.shape)
        E_ms_a_expanded = E_ms_a [:, :, np.newaxis]

        for e_ms_n in E_ms_n_list:
            print (f"ccount = {ccount}", flush=True)
            ccount += 1

            ax = axes[rcount][ccount]
            ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', labelsize=3, pad=2)
            ax.tick_params(axis='x', labelsize=3, width=2)
            ax.tick_params(axis='y', labelsize=3, width=2)


            E_ms_n = e_ms_n * np.ones (E_mm_a.shape)
            E_ms_n_expanded = E_ms_n[:, :, np.newaxis]

            print ("Calculating chis...", flush=True)
            Z_expanded = chi (E_mm_n_expanded, E_mm_a_expanded, E_ms_n_expanded, E_ms_a_expanded, G_expanded, PV_expanded, T_broadcast)
            print ("Calculated!", flush=True)

            print ("Processing the calculated chis...", flush=True)
            hold      = (np.max(Z_expanded, axis=-1)>0) & (Z_expanded[:,:,0]<0)
            Z_cg      = np.where (hold, np.max(Z_expanded, axis=-1) - (Z_expanded[:,:,0]), 0)

            hold   = (np.max(Z_expanded, axis=-1)<0) & (Z_expanded[:,:,0]<0)
            # Z_cc   = np.where (hold, np.max(Z_expanded, axis=-1), 0)    
            Z_cc   = np.where (hold, np.max(Z_expanded, axis=-1) - (Z_expanded[:,:,0]), 0)    

            hold   = (np.min(Z_expanded, axis=-1)>0) & (Z_expanded[:,:,0]>0)
            # Z_gg   = np.where (hold, np.min(Z_expanded, axis=-1) , 0)    
            Z_gg   = np.where (hold, (Z_expanded[:,:,0]) - np.min(Z_expanded, axis=-1) , 0)    

            hold   = (np.min(Z_expanded, axis=-1)<0) & (Z_expanded[:,:,0]>0)
            # Z_gc   = np.where (hold, np.min(Z_expanded, axis=-1) , 0)    
            Z_gc   = np.where (hold, (Z_expanded[:,:,0]) - np.min(Z_expanded, axis=-1) , 0)             

            print ("Processed!", flush=True)

            print ("begin plotting...", flush=True)
            try:
                
                if rcount != ccount:
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_cg, cmap="Reds"    ,   norm=colors.LogNorm(vmin=np.min(Z_cg[Z_cg>0]), vmax=np.max(Z_cg[Z_cg>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_cc, cmap="Greens"  ,   norm=colors.LogNorm(vmin=np.min(Z_cc[Z_cc>0]), vmax=np.max(Z_cc[Z_cc>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_gg, cmap="Blues"   ,   norm=colors.LogNorm(vmin=np.min(Z_gg[Z_gg>0]), vmax=np.max(Z_gg[Z_gg>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_gc, cmap="Purples" ,   norm=colors.LogNorm(vmin=np.min(Z_gc[Z_gc>0]), vmax=np.max(Z_gc[Z_gc>0])), shading="auto" )

                else:
                    # ax.pcolormesh ( E_mm_n, E_mm_a, Z_cg, cmap="Reds"    ,   norm=colors.LogNorm(vmin=np.min(Z_cg[Z_cg>0]), vmax=np.max(Z_cg[Z_cg>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_cc, cmap="Greens"  ,   norm=colors.Normalize(vmin=np.min(Z_cc[Z_cc>0]), vmax=np.max(Z_cc[Z_cc>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_gg, cmap="Blues"   ,   norm=colors.Normalize(vmin=np.min(Z_gg[Z_gg>0]), vmax=np.max(Z_gg[Z_gg>0])), shading="auto" )
                    ax.pcolormesh ( E_mm_n, E_mm_a, Z_gc, cmap="Purples" ,   norm=colors.Normalize(vmin=np.min(Z_gc[Z_gc>0]), vmax=np.max(Z_gc[Z_gc>0])), shading="auto" )
            
            except ValueError:
                pass
                # ax.pcolormesh (E_mm_n, E_mm_a, Z, cmap='Reds', vmin=0, vmax=1)   

            del E_ms_n
            del E_ms_n_expanded
            del Z_expanded
            del Z_cg
            del Z_cc
            del Z_gg
            del Z_gc
            del hold

            ax.set_xlim (-plot_lim, plot_lim)
            ax.set_ylim (-plot_lim, plot_lim)
            ax.set_yticks([-plot_lim,0,plot_lim])
            ax.set_xticks([-plot_lim,0,plot_lim])
            ax.set_xlabel ("$\\epsilon _{mm} ^{n}$", fontsize=3, weight='bold', labelpad=2)
            ax.set_ylabel ("$\\epsilon _{mm} ^{a}$", fontsize=3, labelpad=2)
            ax.set_xticklabels (ax.get_xticks(), weight='bold')
            ax.set_yticklabels (ax.get_yticks(), weight='bold')
            ax.minorticks_on()
            print ("plotted!", flush=True)

        del E_ms_a
        del E_ms_a_expanded

    fig.savefig ("fast-plots-multicolor-medium-v2.png", dpi=1200, bbox_inches="tight")

    stop = time.time()
    print(f"Time required by this computation is {stop-start} seconds.")


