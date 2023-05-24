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

    frac_a = lambda x_s, phi_b: x_s * (1-phi_b)

    Tup    = lambda phi_a, phi_b, E_ab, E_bc, E_ac: 1/ (1+(N-1)*phi_b) * \
    ( E_ac * phi_a - E_ac * phi_a **2 + N * E_bc * phi_b + N* E_ab * phi_a * phi_b - E_ac * phi_a * phi_b - N * E_bc * phi_a * phi_b - N * E_bc* phi_b ** 2 + 
        0.5* np.sqrt ( -4 * N * (E_ab**2 + (E_ac - E_bc)**2 - 2 * E_ab * (E_ac + E_bc) ) * phi_a * phi_b * (phi_a + phi_b - 1) * ( 1 + phi_b * (N-1) ) + \
            4 * (E_ac * phi_a * (phi_a + phi_b - 1) + N * phi_b * (-E_ab * phi_a + E_bc * (-1 + phi_a + phi_b) ) ) ** 2 ) )


    E_ab  = 1        # a is a bad solvent
    E_bc  = 1        # c is a bad solvent
    phi_b = 0.01 
    x_s_range   = np.linspace (0, 1, 100) 

    E_ac_list = [-2, -1, -0.5, 0, 0.5, 1, 2 ]#, 2, 3, 4, 5] # , 6, 7, 8, 10, 25] # change up how things are mixing 

    for idx,E_ac in enumerate(E_ac_list): 
        
        T_specific = lambda x_s: Tup (1-frac_a (x_s, phi_b), phi_b, E_ab, E_bc, E_ac)
        T_solution = T_specific (x_s_range)
        x_s_range  = x_s_range [~np.isnan(T_solution)]
        T_solution = T_solution[~np.isnan(T_solution)]
        p = ax.plot (1-x_s_range, T_solution, ls='-', lw=1, zorder=10, solid_capstyle='round', label=f"$\\chi _{{sc}}  = {E_ac}/T$") 
        if np.max (T_solution) > 0:
            ax.axvline (x=1-x_s_range[np.argmax(T_solution)], ls='--', lw=0.5, color=p[0].get_color())
        p.clear()


    ax.minorticks_on ()
    ax.set_xticks (np.linspace (0, 1, 6))
    ax.legend(fontsize=4, loc="upper right")
    # ax.set_ylim (bottom=0)
    ax.set_xlim (0, 1)
    plt.savefig ("plots-cosolvency-solventmixing.png", dpi=1200) # bbox_inches='tight', dpi=1200)
