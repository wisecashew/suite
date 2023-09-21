import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.colors as colors 
from scipy.optimize import fsolve
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator, AutoMinorLocator, MultipleLocator
from scipy.spatial.distance import cdist
import sys
import argparse
import linecache
np.set_printoptions(threshold=sys.maxsize)
import warnings 


import argparse 

parser = argparse.ArgumentParser(description="Locate the /c/ritical points on the /spin/odal diagram. This program will create one plot, and you can customize what you want on the plot.")
parser.add_argument('-vs',      metavar='vs',      dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',      dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',      dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--no-rtw',              dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--img-name',            dest='img',       action='store', type=str,  default="None", help='name of the image to be created (default: all of the inputs in the imagename).')
args = parser.parse_args()

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    if args.nrtw:
        return f"beep.\n"
    else:
        return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

if __name__=="__main__":

    print ("Revving up the program...", flush=True)

    lsize = 3
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(num=1, figsize=(6,6))
    ax  = fig.add_subplot (projection="3d")

    vs     = args.vs
    vc     = args.vc
    vp     = args.vp


    discriminant = lambda phi_s, chi_ps, chi_pc, chi_sc: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc)) + (vp - 2*phi_s*vp *vs *chi_ps + vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc) ) ) )**2 # -4* vb * (- vc + phi_s * vc - phi_s * vs + 2*phi_s*vc*vs*chi_sc - 2*phi_s**2*vc*vs*chi_sc) * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) **2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) + \
    # (-1 + 2 * phi_s * chi_sc + vb * (1 - 2*chi_pc - phi_s * (chi_ps ** 2 + chi_sc **2 - 2*chi_sc*chi_pc + (chi_pc -2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) + phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) ) ) ** 2
    denom  = lambda phi_s, chi_ps, chi_pc, chi_sc: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))  # 1 / (2*vb * (2*chi_pc + phi_s * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2*chi_ps * (chi_sc + chi_pc) ) ) )
    prefac = lambda phi_s, chi_ps, chi_pc, chi_sc: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp * (2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc) ) ) # 1 - 2 * phi_s * chi_sc + vb * ( -1 + 2 * chi_pc + phi_s * (chi_ps ** 2 + chi_sc ** 2 - 2*chi_sc * chi_pc + (chi_pc - 2) * chi_pc - 2 * chi_ps * (-1 + chi_sc + chi_pc) ) - phi_s ** 2 * (chi_ps ** 2 + (chi_sc - chi_pc) ** 2 - 2 * chi_ps * (chi_sc + chi_pc) ) )
    root_up  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) + np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )
    root_lo  = lambda phi_s, chi_ps, chi_pc, chi_sc: denom (phi_s, chi_ps, chi_pc, chi_sc) * ( prefac (phi_s, chi_ps, chi_pc, chi_sc) - np.sqrt(discriminant (phi_s, chi_ps, chi_pc, chi_sc) ) )

    def stab_crit (p_s, p_p, c_ps, c_pc, c_sc):
        vals = (1/(vp*p_p) + 1/(vc*(1-p_s - p_p)) - 2 * c_pc[:, :, :, np.newaxis]) * (1/(vs*p_s) + 1/(vc*(1-p_s - p_p)) - 2 * c_sc[:, :, :, np.newaxis]) - (1/(vc*(1-p_s-p_p)) + c_ps[:, :, :, np.newaxis] - c_pc[:, :, :, np.newaxis] - c_sc[:, :, :, np.newaxis]) ** 2
        negs = np.any(vals < 0, axis=-1)
        unstab = [c_ps[negs],  c_pc[negs],  c_sc[negs] ]
        stab   = [c_ps[~negs], c_pc[~negs], c_sc[~negs]]
        return unstab, stab

    # get all the volume fractions
    p_s_space = np.arange (0.001, 1-0.001, 0.001)
    p_s = np.repeat (p_s_space, len(p_s_space))

    p_p = np.zeros (p_s.shape)
    for i in range (len(p_s_space)):
        p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

    linrange = 20
    
    sc_range = np.linspace(-15, 0,  20)
    ps_range = np.linspace(-15, 0,  20)
    pc_range = np.linspace(-15, 0,  20)
    chisc_range, chips_range, chipc_range = np.meshgrid(sc_range, ps_range, pc_range)

    c_unstab = "coral"
    c_stab   = "grey"

    unstab, stab = stab_crit (p_s, p_p, chisc_range, chips_range, chips_range)
    print (f"vals.shape = {vals.shape}")

    ax.scatter (unstab[0].reshape(-1), unstab[1].reshape(-1), unstab[2].reshape(-1), marker='o', c=c_unstab)
    ax.scatter (stab[0].reshape(-1), stab[1].reshape(-1), stab[2].reshape(-1), marker='o', c=c_stab, alpha=0.2)

    ax.set_xlim (-15, 0)
    ax.set_ylim (-15, 0)
    ax.set_zlim (-15, 0)

    ax.set_xlabel("$\\chi _{ps}$")
    ax.set_ylabel("$\\chi _{pc}$")
    ax.set_zlabel("$\\chi _{sc}$")

    ax.grid()
    plt.show()


    if args.img != "None":
        if "." in args.img:
            plt.savefig (args.img+".png", dpi=1200)
        else:
            plt.savefig (args.img, dpi=1200)
    else:
        plt.savefig (f"3d_mapping.png", dpi=1200)
    
    print ("Completed heat map computation.")

