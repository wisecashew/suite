#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import matplotlib.patheffects as pe 
from scipy.optimize import fsolve
from scipy.optimize import root
from matplotlib.ticker import StrMethodFormatter
import scipy.optimize as opt 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import Locator


if __name__=="__main__":

    EDIT = True

    lsize = 8
    plt.rcParams['font.family'] = 'Arial'
    font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    fig = plt.figure(figsize=(2.5,2.0), constrained_layout=True)
    fig.tight_layout()
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=lsize, pad=5)
    ax.tick_params(axis='y', labelsize=lsize)

    
    chi = lambda E, T: E/T
    N   = 100

    frac_a = lambda r, phi_b: (1-phi_b)*r/(r+1)

    Tup    = lambda phi_a, phi_b, E_ab, E_bc, E_ac: 1/ (1+(N-1)*phi_b) * \
    ( E_ac * phi_a - E_ac * phi_a **2 + N * E_bc * phi_b + N* E_ab * phi_a * phi_b - E_ac * phi_a * phi_b - N * E_bc * phi_a * phi_b - N * E_bc* phi_b ** 2 + 
        0.5* np.sqrt ( -4 * N * (E_ab**2 + (E_ac - E_bc)**2 - 2 * E_ab * (E_ac + E_bc) ) * phi_a * phi_b * (phi_a + phi_b - 1) * ( 1 + phi_b * (N-1) ) + \
            4 * (E_ac * phi_a * (phi_a + phi_b - 1) + N * phi_b * (-E_ab * phi_a + E_bc * (-1 + phi_a + phi_b) ) ) ** 2 ) )

    E_ab  = -150
    E_bc  = -500
    E_ac  = 0
    phi_b_range = np.arange (0.01, 0.5, 0.005)

    ratio_list = [1, 0.5, 2, 0.2, 5, 0.1, 10, 0.02, 50]

    # color_dict ={1: "steelblue", 2: "coral", 0.5: "coral", 5: "pink", 0.2: "pink", 0.1: "forestgreen", 10: "forestgreen", 0.02: "darkred", 50: "darkred"}

    for ratio in ratio_list: 
        
        T_specific = lambda phi_b: Tup (frac_a (ratio, phi_b), phi_b, E_ab, E_bc, E_ac) 
        T_solution = T_specific (phi_b_range)     
        # if ratio <= 1:
        ax.plot (phi_b_range, T_solution, ls='-', lw=1, zorder=10, solid_capstyle='round', label=f"r={ratio}") # , label=f"r = {ratio}, {1/ratio}") , c=color_dict[ratio])

        # else:
        #     ax.plot (phi_b_range, T_solution, ls='-', lw=1, zorder=10, solid_capstyle='round') # , label="_nolabel_", c=color_dict[ratio])            
    


    ax.minorticks_on ()
    ax.set_xticks (np.linspace (0, 0.5, 6))
    ax.legend(fontsize=4, loc="upper right")
    # ax.set_ylim (bottom=0, top=500)
    # ax.set_xlim (0, 0.5)
    if EDIT:
        plt.savefig ("plots-cononsolvency-fig4-edit.png", dpi=1200) # bbox_inches='tight', dpi=1200)
    else:
        plt.savefig ("plots-cononsolvency-fig4.png", dpi=1200) # bbox_inches='tight', dpi=1200)
