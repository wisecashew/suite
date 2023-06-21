import numpy as np
import matplotlib 
import matplotlib.font_manager
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.patheffects as pe 
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
import argparse 
from pathlib import Path


parser = argparse.ArgumentParser (description="Plots phase diagrams.")
parser.add_argument ("--png-name", dest="pn", type=str, action='store', help="Name of image to be made.")
parser.add_argument ("--pv", dest="pv", type=np.float128, action='store', help="Value of pv.")

args = parser.parse_args() 


if __name__=="__main__":

    fpath = Path(matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf") 
    lsize = 8
    fdict = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(figsize=(2.5,2.0), constrained_layout=True)
    # fig.tight_layout()
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=lsize)
    ax.tick_params(axis='y', labelsize=lsize)

    g    = 0.25
    zmm  = lambda emma, emmn, T: g*np.exp (-1/T * emma, dtype=np.float128) + (1-g)*np.exp (-1/T * emmn, dtype=np.float128)
    zms  = lambda emsa, emsn, T: g*np.exp (-1/T * emsa, dtype=np.float128) + (1-g)*np.exp (-1/T * emsn, dtype=np.float128)
    fmma = lambda emma, emmn, T: g*np.exp (-1/T * emma, dtype=np.float128)/zmm(emma, emmn, T)
    fmsa = lambda emsa, emsn, T: g*np.exp (-1/T * emsa, dtype=np.float128)/zms(emsa, emsn, T)

    def chi (emma, emmn, emsa, emsn, pv, T):
        t1 = pv*(fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn) + (1-pv)*emmn
        return 24/T * (t1 - 0.5 * t2)

    def rhs (phi, _m):
        return 1/2 * (1/(1-phi) + 1/(_m*phi))

    hexcolor_cg = '#B91F72'
    hexcolor_cc = '#369DE8'
    hexcolor_gg = '#1FB967' 
    hexcolor_gc = '#B9B41F'
    cols = [hexcolor_cg, hexcolor_gg, hexcolor_cc, hexcolor_gc]
    # blue, red, green -> green, blue, red
    T         = np.logspace (-2, 2, 200)
    colours   = ["#33CC37", "#3733CC", "#CC3733", "#936C77", "#936C77", "#936C77"]
    emsa_list = [0] # [-2, -1.5, -1, -0.5, -0.1, 0]
    elow    = np.min(emsa_list)
    ehigh   = np.max (emsa_list)
    pv      = args.pv
    ecenter = (elow+ehigh)/2

    norm = matplotlib.colors.Normalize (vmin=-4,  vmax=0)

    eparams = [(-1,-1,-0.8,0), (-1,-1,0,0), (-1,-1,-1,-1), (-2.2,0,-1,-1)]

    for idx,pars in enumerate(eparams):
        # rgba_color   = cm.rainbow( norm (pars) )
        # ax.plot (T, chi(_emma, _emmn, _emsa, _emsn, pv, T), ls='-', lw=1, c=rgba_color, zorder=10, solid_capstyle='round')
        ax.plot (T, chi(pars[0], pars[1], pars[2], pars[3], pv, T), ls='-', lw=1, c=cols[idx], zorder=10, solid_capstyle='round')

    # ax.minorticks_on ()
    # ax.set_xticks (np.logspace (-2, 2, 5))
    ax.set_xscale ("log")
    ax.set_yscale ("symlog")
    # ax.set_yticks ([0.1, 1, 10, 50 ])
    ax.set_yticklabels (ax.get_yticks(), font=fpath, fontdict=fdict)
    if pv == 1.0:
        ax.set_xticklabels (ax.get_xticks(), fontdict=fdict, font=fpath)
        plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    else:
        ax.set_xticklabels([])

    ax.axhline (y=0, c='k', linestyle='--', lw=0.5)

    ax.set_ylim (-1000, 1000)
    ax.set_xlim (0.01, 100)
    plt.savefig (args.pn+".png", dpi=1200) # bbox_inches='tight', dpi=1200)

