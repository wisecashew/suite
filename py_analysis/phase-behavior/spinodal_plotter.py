import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mplc
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import fsolve
from matplotlib.ticker import StrMethodFormatter
import time
import warnings
import linecache
import argparse
from pathlib import Path


parser = argparse.ArgumentParser (description="Plot the spinodals given a certain set of energetic parameters. You can also input the transition temperature obtained via simulation to see where it lies wrt to the spinodal.")
parser.add_argument("--EMMA", dest="emma", type=float, action='store', help="Enter monomer-monomer aligned interactions (default: 0).", default=0)
parser.add_argument("--EMMN", dest="emmn", type=float, action='store', help="Enter monomer-monomer misaligned interactions (default: 0).", default=0)
parser.add_argument("--EMSA", dest="emsa", type=float, action='store', help="Enter monomer-solvent aligned interactions (default: 0).", default=0)
parser.add_argument("--EMSN", dest="emsn", type=float, action='store', help="Enter monomer-solvent misaligned interactions (default: 0).", default=0)
parser.add_argument("--ESSA", dest="essa", type=float, action='store', help="Enter solvent-solvent aligned interactions (default: 0).", default=0)
parser.add_argument("--ESSN", dest="essn", type=float, action='store', help="Enter solvent-solvent misaligned interactions (default: 0).", default=0)
parser.add_argument("--PV", dest="pv", type=float, action='store', help="Enter the operational volume (default: 1).", default=1)
parser.add_argument("--PW", dest="pw", type=float, action='store', help="Enter the microcanonical probability (default: 0.25).", default=0.25)
parser.add_argument("--dop",  dest="dop",  type=float, action='store', help="Degree of polymerization of the polymer.")
parser.add_argument("--T", dest="T", type=float, nargs='+', action='store', help="Temperatures at which we see some kind of transition in the single chain simulation.")
parser.add_argument("--hexcode", dest="hex", type=str, action='store', help="Hex code used to color the spinodal. You need to enclose the hexcode in double quotation marks when using it as in input e.g. `--hexcode \"#1AC6E5\"`.")
parser.add_argument("--png-name", dest="pn", type=str, action='store', help="Name of image to be made. If you have a period (`.`), the program will add a .png extension.")

args = parser.parse_args() 

def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"Probably a math sqrt error (line {lineno}).\n"

warnings.formatwarning = custom_warning_format


# Define the color gradients and their corresponding positions
c1 = "#1AC6E5"
c2 = "#244BDB"
c3 = "#C157DB"
c4 = "#D95B59"
c5 = "#BFDA58"
c6 = "#2ED16A"

colors = [c1, c2, c2, c3, c3, c4, c4, c5, c5, c6]
positions_raw   = np.array([-1, -0.53, -0.53, -0.525, -0.5250, -0.5235, -0.5235, -0.5, -0.5, 0])
k = 1/(np.max(positions_raw) - np.min (positions_raw))
d = 0-np.min(positions_raw)
positions_clean = k * ( positions_raw + d )

# Create the colormap
cmap = cm.rainbow # LinearSegmentedColormap.from_list('my_colormap', list(zip(positions_clean, colors)))

norm = mplc.Normalize (vmin=-4, vmax=0)


def good_transitions (arr):

    positive = arr[0] > 0

    if positive:
        pass
    else:
        return False

    transition = 0
    for i in range(1, len(arr)):
        current_positive = arr[i] > 0

        if current_positive != positive:
            transition += 1

        positive = current_positive
        
    if transition >= 2:
        return True
    else:
        return False


zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float128)
zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float128)
zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128) + (1-pw)*np.exp (-1/T * essn, dtype=np.float128)
fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128)/zmm(emma, emmn, pw, T)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128)/zms(emsa, emsn, pw, T)
fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128)/zss(essa, essn, pw, T)

###############################################################

class Phase:

    def __init__ (self, param_list, N):
        self.EMSA   = param_list[0]
        self.EMSN   = param_list[1]
        self.EMMA   = param_list[2]
        self.EMMN   = param_list[3]
        self.ESSA   = param_list[4]
        self.ESSN   = param_list[5]
        self.PV     = param_list[6]
        self.PWMS   = param_list[7]
        self.PWMM   = param_list[8]
        self.PWSS   = param_list[9]
        self.N      = N
        return

    def reset_params (self, param_list):
        self.EMSA   = param_list[0]
        self.EMSN   = param_list[1]
        self.EMMA   = param_list[2]
        self.EMMN   = param_list[3]
        self.ESSA   = param_list[4]
        self.ESSN   = param_list[5]
        self.PV     = param_list[6]
        self.PWMS   = param_list[7]
        self.PWMM   = param_list[8]
        self.PWSS   = param_list[9]
        return

    def print_params (self):
        print (f"The parameters of this curve are EMSA = {self.EMSA}, EMSN = {self.EMSN}, EMMA = {self.EMMA}, EMMN = {self.EMMN}, ESSA = {self.ESSA}, ESSN = {self.ESSN}, \
PV = {self.PV}, PWMS = {self.PWMS}, PWMM = {self.PWMM}, PWSS = {self.PWSS}", flush=True)
        return

    def chi (self, T):
        c = 24 * (self.PV * ( (fmsa (self.EMSA, self.EMSN, self.PWMS, T) * self.EMSA + (1 - fmsa (self.EMSA, self.EMSN, self.PWMS, T) ) * self.EMSN) - 1/2 * \
        ( (fmma (self.EMMA, self.EMMN, self.PWMM, T) * self.EMMA + (1-fmma (self.EMMA, self.EMMN, self.PWMM, T) ) * self.EMMN) + \
          (fssa (self.ESSA, self.ESSN, self.PWSS, T) * self.ESSA + (1-fssa (self.ESSA, self.ESSN, self.PWSS, T) ) * self.ESSN) ) )
        + (1-self.PV) * (self.EMSN - 1/2 * (self.EMMN + self.ESSN) ) ) / T

        return np.array(c, dtype=np.float128)

    def spinodal (self, T):
        p1 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) - np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        p2 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) + np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        return (np.array(p1, dtype=np.float128), np.array(p2, dtype=np.float128), T)

    def delta_phi (self, T):
        p1 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) - np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        p2 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) + np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        return p1 - p2


    # end of class Phase 

###############################################################

if __name__=="__main__":

    start = time.time()
    fpath = Path (matplotlib.get_data_path(), "/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis/arial.ttf")
    lsize = 12
    fdict = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

    energy_param_list  = [args.emsa, args.emsn, args.emma, args.emmn, args.essa, args.essn, args.pv, args.pw, args.pw, args.pw]
    N                  = args.dop
    PhaseDiag          = Phase (energy_param_list, N)
    T                  = np.hstack((np.logspace (-3, np.log10(30), int(1e+7) ), np.linspace (0.05,0.15, int(1e+6)), np.linspace (0.1, 1.0, int(1e+6) ) ) )
    T                  = np.sort (T, kind="mergesort")


    fig = plt.figure(figsize=(2.5,2.5), constrained_layout=True)
    fig.tight_layout()
    ax  = plt.axes ()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=lsize, pad=5)
    ax.tick_params(axis='y', labelsize=lsize)

    PhaseDiag.print_params()
    if args.hex == None:
        rgba_color   = cmap(norm(args.emsa))
    else:
        rgba_color = args.hex

    spinodal = PhaseDiag.spinodal (T)
    T1       = spinodal[2][spinodal[0] < 1]
    arm1     = spinodal[0][spinodal[0] < 1]
    T1       = T1  [arm1 > 0]
    arm1     = arm1[arm1 > 0]
    line = ax.plot(arm1, T1, lw=1.0, markersize=0, c=rgba_color, solid_capstyle="round",label="_nolabel_")
    T1       = spinodal[2][spinodal[1] < 1]
    arm1     = spinodal[1][spinodal[1] < 1]
    T1       = T1  [arm1 > 0]
    arm1     = arm1[arm1 > 0]
    ax.plot  (arm1, T1, lw=1.0, markersize=0, color=rgba_color, solid_capstyle="round", label="analytic")

    T_sim = args.T
    if T_sim == None:
        print (f"No transition temperatures provided.")
    else:
        for idx,T_tr in enumerate(T_sim):
            if idx == 0:
                ax.axhline (y=T_tr, c='k', ls='--', label="from simulation")
            else:
                ax.axhline (y=T_tr, c='k', ls='--', label="_nolabel_")


    del T1
    del arm1
    del spinodal

    ax.legend (loc="upper right", fontsize=4.5)
    ax.set_yscale ("log")
    ax.set_ylim (1e-3, 50)
    ax.set_xlim (0, 1)

    ax.set_yticks ([0.001, 0.01, 0.1, 1.0, 10, 50])
    ax.set_xticks ([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels ([0.001, 0.01, 0.1, 1.0, 10, 50], fontdict=fdict, font=fpath)
    ax.set_xticklabels ([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontdict=fdict, font=fpath)

    # make sure the name of the image file is correct
    img_name = args.pn
    if "." in img_name:
        img_name = img_name + ".png"
    else:
        pass

    fig.savefig (img_name, dpi=1200, bbox_inches="tight")
    print (f"Plotted!", flush=True)
    stop = time.time()
    print (f"Elapsed time is {stop-start} seconds.", flush=True)
