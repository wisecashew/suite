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


    g    = 0.5
    zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float128)
    zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float128)
    zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128) + (1-pw)*np.exp (-1/T * essn, dtype=np.float128)
    fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128)/zmm(emma, emmn, pw, T)
    fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128)/zms(emsa, emsn, pw, T)
    fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128)/zss(essa, essn, pw, T)

    def chi (emma, emmn, emsa, emsn, essa, essn, pw, pv, T):
        t1 = pv*(fmsa(emsa, emsn, pw, T)*emsa + (1-fmsa(emsa, emsn, pw, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, pw, T)*emma + (1-fmma(emma, emmn, pw, T))*emmn) + (1-pv)*emmn
        t3 = pv*(fssa(essa, essn, pw, T)*essa + (1-fssa(essa, essn, pw, T))*essn) + (1-pv)*essn
        return 24/T * (t1 - 0.5 * (t2 + t3))

    def rhs (phi, _m):
        return 1/2 * (1/(1-phi) + 1/(_m*phi))

    hexcolor_cg = '#B91F72'
    hexcolor_cc = '#369DE8'
    hexcolor_gg = '#1FB967' 
    hexcolor_gc = '#B91F72' # '#B9B41F'
    cols = [hexcolor_cc, hexcolor_cg, hexcolor_gg, hexcolor_gc]
    # blue, red, green -> green, blue, red
    T         = np.logspace (-5, np.log10(500), 300)
    emsa_list = [0]
    elow      = np.min(emsa_list)
    ehigh     = np.max (emsa_list)
    pv        = args.pv
    ecenter   = (elow+ehigh)/2

    norm = matplotlib.colors.Normalize (vmin=-4,  vmax=0)
    
    eparams = [(0,0,-1,-1,0,0), (-1,-1,-0.6,0,0,0)]# [(-2.001, -2.001, -1, 0, 0, 0)] # , (-1, -1, 0, 0, 0, 0), (0, 0, -1, -1, 0, 0), (-1, 0, -0.3, -0.3, 0, 0) ]
    eparams = [(-0.01,-0.01,0,0,0,0),(-2.001/10,-2.001/10,-1/10,0,0,0)]# [(-2.001, -2.001, -1, 0, 0, 0)] # , (-1, -1, 0, 0, 0, 0), (0, 0, -1, -1, 0, 0), (-1, 0, -0.3, -0.3, 0, 0) ]

    fig = plt.figure(num=0, figsize=(2.5,2.0), constrained_layout=True)
    ax  = plt.axes ()
    ax.axis ('off')
    ax.tick_params(direction='in', bottom=False, top=False, left=False, right=False, which='both')
    for idx,pars in enumerate(eparams):
        Y = chi(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], g, pv, T)
        ax.plot (T, chi(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], g, pv, T), ls='-', lw=1, c=cols[idx+2], zorder=10, solid_capstyle='round')

        ax.set_xscale ("log")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axhline (y=0, c='k', linestyle='-', lw=1, zorder=11)
        ax.axhline (y=0.6924016952966369, c='r', linestyle='--', lw=1, zorder=11)
        ax.axvline (x=0, c='k', linestyle='-',  lw=1, zorder=11)
        ax.set_ylim (-3, 10)
        ax.set_xlim (0.0001, 500)


    plt.savefig (args.pn, dpi=1200) # bbox_inches='tight', dpi=1200)

