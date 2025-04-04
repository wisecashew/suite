#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
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
parser.add_argument ("--pv", dest="pv", type=np.float64, action='store', help="Value of pv.")

args = parser.parse_args() 


if __name__=="__main__":

    fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
    lsize = 10
    fdict = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(figsize=(2.5,2.0), constrained_layout=True)
    fig.tight_layout()
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=0, pad=5)
    ax.tick_params(axis='y', labelsize=lsize)

    g    = 0.25
    zmm  = lambda emma, emmn, T: g*np.exp (-1/T * emma, dtype=np.float64) + (1-g)*np.exp (-1/T * emmn, dtype=np.float64)
    zms  = lambda emsa, emsn, T: g*np.exp (-1/T * emsa, dtype=np.float64) + (1-g)*np.exp (-1/T * emsn, dtype=np.float64)
    fmma = lambda emma, emmn, T: g*np.exp (-1/T * emma, dtype=np.float64)/zmm(emma, emmn, T)
    fmsa = lambda emsa, emsn, T: g*np.exp (-1/T * emsa, dtype=np.float64)/zms(emsa, emsn, T)

    def chi (emma, emmn, emsa, emsn, pv, T):
        t1 = pv*(fmsa(emsa, emsn, T)*emsa + (1-fmsa(emsa, emsn, T))*emsn) + (1-pv)*emsn
        t2 = pv*(fmma(emma, emmn, T)*emma + (1-fmma(emma, emmn, T))*emmn) + (1-pv)*emmn
        return 24/T * (t1 - 0.5 * t2)

    def rhs (phi, _m):
        return 1/2 * (1/(1-phi) + 1/(_m*phi))



    # parameter list
    # phi_list   = np.arange (0.01, 1.0, 0.01)
    phi_list   = np.arange (0.01, 1.0, 0.00001) # np.arange (0.3, 0.36, 0.001) # np.arange (0.01, 1.0, 0.0001)
    seeds      = [0.01, 0.05, 0.1, 1.0, 3.4721485, 5.25, 10.0, 25.0, 50.0, 100.0]
    roots      = np.array([])
    root_error = np.array([])
    T_lower    = []
    T_upper    = []

    # blue, red, green -> green, blue, red
    colours   = ["#33CC37", "#3733CC", "#CC3733", "#936C77", "#936C77", "#936C77"]
    emsa_list = [-4.0, -3.5, -3.25, -3.1, -3, -2.5, -2, -1.75, -1.5, -1.3, -1.235, -1.22, -1.2, -1.1, -1, -0.8, -0.75, -0.65,  -0.5, -0.25, -0.1]
    elow    = np.min(emsa_list)
    ehigh   = np.max (emsa_list)
    pv      = args.pv
    ecenter = (elow+ehigh)/2
    # norm = matplotlib.colors.TwoSlopeNorm (vmin=elow, vcenter=ecenter, vmax=ehigh)
    norm = matplotlib.colors.Normalize (vmin=-4,  vmax=0)


    # energies 
    _emma = -1; _emmn = -1; _emsn = 0; 
    my_phi = []

    for idx, _emsa in enumerate(emsa_list):
        print ("emsa = ", _emsa)
        T_lower.clear ()
        T_upper.clear ()
        my_phi.clear ()

        for _phi in phi_list:

            roots        = np.array([])
            root_error   = np.array([])
            specific_chi = lambda T: chi (_emma, _emmn, _emsa, _emsn, pv, T) - rhs (_phi, 32)

            for s in seeds:
                roots = np.hstack ( ( roots, fsolve (specific_chi, s) ) )
                root_error = np.hstack ( ( root_error, abs (specific_chi(roots[-1]) ) ) )

            hold = roots > 0
            roots = roots [hold]
            root_error = root_error [hold]


            hold = root_error < 1e-9
            roots = roots [hold]
            root_error = root_error [hold]

            # print ("phi_b = ", _phi)
            # print ("roots = ", roots)
            # print (f"np.max(roots) = {np.max(roots)}")
            # print ("root_error = ", root_error)

            try:
                T_lower.append (np.min (roots))
                T_upper.append (np.max (roots))
                my_phi.append  (_phi)

            except ValueError:
                print ("Root finding was unstable.")
                # print ("phi_b = ", _phi)


        rgba_color   = cm.rainbow( norm (_emsa) )
        ax.plot (my_phi[::1], T_lower[::1], ls='-', lw=1, c=rgba_color, zorder=10, solid_capstyle='round', label="$\\epsilon _{ms} ^{\\parallel}$ = " + str(_emsa) ) # , path_effects=[pe.Stroke(linewidth=3.5, foreground='k'), pe.Normal()]) # , clip_on=False)
        ax.plot (my_phi[::1], T_upper[::1], ls='-', lw=1, c=rgba_color, zorder=10, solid_capstyle='round')# , path_effects=[pe.Stroke(linewidth=3.5, foreground='k'), pe.Normal()]) #, clip_on=False)

    ax.minorticks_on ()
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_yscale ("log")
    ax.set_yticks ([0.1, 1, 10, 50 ])
    ax.set_yticklabels (ax.get_yticks(), fontdict=fdict, font=fpath)
    if pv == 1.0:
        ax.set_xticklabels (ax.get_xticks(), fontdict=fdict, font=fpath)
        plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    else:
        ax.set_xticklabels([])

    ax.set_ylim (0.1, 50)
    ax.set_xlim (0, 1)
    plt.savefig (args.pn+".png", dpi=1200)


