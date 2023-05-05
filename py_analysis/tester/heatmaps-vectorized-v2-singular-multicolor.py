#!/home/satyend/.conda/envs/data_analysis/bin/python

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


def map_maker (hex_code):

    rgbcolor    = colors.hex2color (hex_code)
    cmap_colors = [(1,1,1), rgbcolor]
    my_cmap     = colors.LinearSegmentedColormap.from_list ('custom_cmap', cmap_colors)

    return my_cmap 


if __name__=="__main__":

    #41CA27 (green), '#D8CA27' (ochre yellow), '#EE9EFE' (pink), '#00A8FF' (coolblue)
    hexcolor_cg = '#FF8D92' # '#B91F72'
    cmap_cg     = map_maker (hexcolor_cg)

    hexcolor_cc = '#B1FAFF' # '#369DE8'
    cmap_cc     = map_maker (hexcolor_cc)

    hexcolor_gg = '#AEFFB3' # '#1FB967' 
    cmap_gg     = map_maker (hexcolor_gg) 

    hexcolor_gc = '#FFFCA5' # '#B9B41F'
    cmap_gc     = map_maker (hexcolor_gc)

    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black','weight': 'normal', 'size': 14}

    E_ms_a_list = [0]
    E_ms_n_list = [-0.5]

    fig = plt.figure( figsize=(3,3), constrained_layout=True)
    ax  = plt.axes()
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
    g = 0.5
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

            # ax = axes[2-rcount][ccount]
            ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', labelsize=10, pad=5)
            ax.tick_params(axis='x', labelsize=10)
            ax.tick_params(axis='y', labelsize=10)


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
                dchi_min = np.min([np.min (Z_cg[Z_cg>0]), np.min (Z_cc[Z_cc>0]), np.min (Z_gg[Z_gg>0]), np.min (Z_gc[Z_gc>0])])
                dchi_max = np.max([np.max (Z_cg[Z_cg>0]), np.max (Z_cc[Z_cc>0]), np.max (Z_gg[Z_gg>0]), np.max (Z_gc[Z_gc>0])])
            except ValueError:
                dchi_min = np.min([np.min (Z_cc[Z_cc>0]), np.min (Z_gg[Z_gg>0]), np.min (Z_gc[Z_gc>0])])
                dchi_max = np.max([np.max (Z_cc[Z_cc>0]), np.max (Z_gg[Z_gg>0]), np.max (Z_gc[Z_gc>0])])
            
            print (f"dchi_min = {dchi_min}", flush=True)
            print (f"dchi_max = {dchi_max}", flush=True)

            try:
                ax.pcolormesh ( E_mm_n, E_mm_a, Z_cg, cmap=cmap_cg,   norm=colors.LogNorm(vmin=0.1, vmax=6000), shading="auto" )
                ax.pcolormesh ( E_mm_n, E_mm_a, Z_cc, cmap=cmap_cc,   norm=colors.LogNorm(vmin=0.1, vmax=6000), shading="auto" )
                ax.pcolormesh ( E_mm_n, E_mm_a, Z_gg, cmap=cmap_gg,   norm=colors.LogNorm(vmin=0.1, vmax=6000), shading="auto" )
                ax.pcolormesh ( E_mm_n, E_mm_a, Z_gc, cmap=cmap_gc,   norm=colors.LogNorm(vmin=0.1, vmax=6000), shading="auto" )

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
            # ax.set_xlabel ("$\\mathbf{ \\epsilon _{mm} } ^{\\perp}$ ", fontsize=6, weight='bold', labelpad=2)
            # ax.set_ylabel ("$\\mathbf{ \\epsilon _{mm} ^{\\parallel} }$", fontsize=6, labelpad=2)
            ax.set_xticklabels (ax.get_xticks(), fontdict=font)
            ax.set_yticklabels (ax.get_yticks(), fontdict=font)
            ax.minorticks_on()
            print ("plotted!", flush=True)

        del E_ms_a
        del E_ms_a_expanded

    fig.savefig ("fast-plots-real-multicolor-singular-v2.png", dpi=1200, bbox_inches="tight")

    stop = time.time()
    print(f"Time required by this computation is {stop-start} seconds.")


