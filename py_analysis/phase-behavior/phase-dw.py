#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
import mpmath as mp 
import argparse 

parser = argparse.ArgumentParser (description="Plots phase diagrams.")
parser.add_argument ("--png-name", dest="pn", type=str, action='store', help="Name of image to be made.")
parser.add_argument ("--pv", dest="pv", type=np.float128, action='store', help="Value of pv.")

args = parser.parse_args() 

# f           = lambda phi, m, x: phi/m * np.log (phi) + (1-phi)*np.log (1-phi) + x * phi * (1-phi)
# dfdphi      = lambda phi, m, x: -1 + 1/m + (1 - phi)*x - phi*x - np.log (1-phi) + np.log(phi)/m
# d2fdphi2    = lambda phi, m, x: 1/(1-phi) + 1/(m*phi) - 2*x
# d3fdphi3    = lambda phi, m, x: 1/(1-phi)**2 - 1/(m*phi*phi)


if __name__=="__main__":

    fig = plt.figure(figsize=(4/1.6,3/1.6), constrained_layout=True)
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both', pad=5, labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

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
    phi_list   = np.arange (0.01, 1.0, 0.01)
    seeds      = [0.005, 0.01, 0.05, 0.1, 1.0, 10.0, 25.0, 50.0]
    roots      = np.array([])
    root_error = np.array([])
    T_lower    = []
    T_upper    = []

    elow    = -0.6
    ehigh   = -0.4
    pv      = args.pv
    ecenter = (elow+ehigh)/2

    emsa_list = [-0.6, -0.4] # np.arange (elow, ehigh+0.05, 0.05)
    norm = matplotlib.colors.TwoSlopeNorm (vmin=elow, vcenter=ecenter, vmax=ehigh)


    # energies 
    _emma = -1; _emmn = -1; _emsn = 0; 
    my_phi = []

    for _emsa in emsa_list:
        print ("emsa = ", _emsa)
        T_lower.clear ()
        T_upper.clear ()
        my_phi.clear ()

        for _phi in phi_list:

            roots      = np.array([])
            root_error = np.array([])
            specific_chi = lambda T: chi (_emma, _emmn, _emsa, _emsn, pv, T) - rhs (_phi, 32)
            

            for s in seeds:
                roots = np.hstack ( ( roots, fsolve (specific_chi, s) ) )
                root_error = np.hstack ( ( root_error, abs (specific_chi(roots[-1]) ) ) )
            
            hold = roots > 0
            roots = roots [hold]
            root_error = root_error [hold]

            hold = root_error < 1e-6
            roots = roots [hold]
            root_error = root_error [hold]

            # print ("phi = ", _phi)    
            # print ("roots = ", roots)
            # print ("root_error = ", root_error)
            
            try:
                T_lower.append (np.min (roots))
                T_upper.append (np.max (roots))
                my_phi.append  (_phi)
            except ValueError:
                print ("Root finding was unstable.")

            
        rgba_color   = cm.PiYG( norm (_emsa) )
        # if _emsa > -1.45:
        #     my_phi.insert  (0, 0)
        #     T_lower.insert (0, 0)
        #     T_upper.insert (0, 0)
        plt.plot (my_phi[::3], T_lower[::3], marker='o', markeredgecolor='k', c=rgba_color)
        plt.plot (my_phi[::3], T_upper[::3], marker='o', markeredgecolor='k', c=rgba_color)

    ax.minorticks_on ()
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.set_yscale ("log")
    ax.set_yticks ([0.1, 0.5, 1, 5, 10, 50])
    ax.set_yticklabels (ax.get_yticks(), weight='bold')
    ax.set_xticklabels (ax.get_xticks(), weight='bold')
    ax.set_yticklabels (["$\\mathbf{0.1}$", "$\\mathbf{0.5}$", "$\\mathbf{1.0}$", "$\\mathbf{5.0}$", "$\\mathbf{10}$", "$\\mathbf{50}$"], weight='bold')
    # plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:1.1f}'))
   
    plt.savefig (args.pn+".png", bbox_inches='tight', dpi=1200)



