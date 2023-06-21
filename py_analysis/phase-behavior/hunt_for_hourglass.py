import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import time

zmm  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float64) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float64)
zms  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float64) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float64)
zss  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float64) + (1-pw)*np.exp (-1/T * essn, dtype=np.float64)
fmma = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float64)/zmm(emma, emmn, pw, T)
fmsa = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float64)/zms(emsa, emsn, pw, T)
fssa = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float64)/zss(essa, essn, pw, T)


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
        self.N    = N
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
        print (f"EMSA = {self.EMSA}, EMSN = {self.EMSN}, EMMA = {self.EMMA}, EMMN = {self.EMMN}, ESSA = {self.ESSA}, ESSN = {self.ESSN}, \
            PV = {self.PV}, PWMS = {self.PWMS}, PWMM = {self.PWMM}, PWSS = {self.PWSS}")
        return

    def chi (self, T):
        c = 24 * (self.PV * ( (fmsa (self.EMSA, self.EMSN, self.PWMS, T) * self.EMSA + (1 - fmsa (self.EMSA, self.EMSN, self.PWMS, T) ) * self.EMSN) - 1/2 * \
        ( (fmma (self.EMMA, self.EMMN, self.PWMM, T) * self.EMMA + (1-fmma (self.EMMA, self.EMMN, self.PWMM, T) ) * self.EMMN) + \
          (fssa (self.ESSA, self.ESSN, self.PWSS, T) * self.ESSA + (1-fssa (self.ESSA, self.ESSN, self.PWSS, T) ) * self.ESSN) ) )
        + (1-self.PV) * (self.EMSN - 1/2 * (self.EMMN + self.ESSN) ) ) / T

        return np.array(c)

    def spinodal (self, T):
        p1 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) - np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        p2 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) + np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        return (np.array(p1), np.array(p2), T)

    def delta_phi (self, T):
        p1 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) - np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        p2 = -1/(4 * self.N * self.chi (T)) * (-1 + self.N - 2 * self.N * self.chi (T) + np.sqrt (-8 * self.N * self.chi (T) + (1 - self.N + 2 * self.N * self.chi (T) ) ** 2 ) )
        return p2 - p1


    # end of class Phase 


if __name__=="__main__":

    start = time.time()
    param_list  = [-5.4, -0, -3, -0, -8.5, 0, 1, 0.01, 0.01, 0.01]
    N           = 100
    PhaseDiag   = Phase (param_list, N)
    PhaseDiag.print_params()
    T           = np.logspace (-2, 2, int(1e+5) )

    # fig0 = plt.figure ()
    # ax0  = plt.axes   ()
    # phi  = np.logspace (-3,np.log10(0.9), int(1e+5))
    # spin_phi = lambda p: 1/(N*p) + 1/(1-p)
    # ax0.axhline (y=np.min(spin_phi(phi)), ls='--')

    # ax0.plot (phi, spin_phi(phi), markersize=0, c="teal")
    # yticks = list (ax0.get_yticks())
    # yticks.append (np.min(spin_phi(phi) ) )
    # ax0.set_yticks (yticks)
    # fig0.savefig ("spin_phi", dpi=1200, bbox_inches="tight")

    fig1 = plt.figure ()
    ax1  = plt.axes   ()
    ax1.set_xlim   (np.min(T), np.max(T))
    ax1.set_xscale ("log")

    chi  = PhaseDiag.chi(T)
    Tchi = T  [~np.isnan(chi)] 
    chi  = chi[~np.isnan(chi)]

    print (f"chi_min = {np.min(chi)}")

    print (f"Creating chi plots...")
    ax1.plot (Tchi, chi, lw=1, markersize=0, c="steelblue")
    fig1.savefig ("chiplots", dpi=1200, bbox_inches="tight")
    print (f"Done!")

    fig2 = plt.figure ()
    ax2  = plt.axes   ()

    print (f"Finding spinodals...")
    spinodal = PhaseDiag.spinodal (T)
    print (f"Found!")

    T1       = spinodal[2][spinodal[0] < 1]
    arm1     = spinodal[0][spinodal[0] < 1]
    T1       = T1   [arm1 > 0]
    arm1     = arm1 [arm1 > 0]
    
    T2       = spinodal[2][spinodal[1] < 1]
    arm2     = spinodal[1][spinodal[1] < 1]
    T2       = T2   [arm2 > 0]
    arm2     = arm2 [arm2 > 0]

    print (f"Plotting spinodal...")
    # ax2.set_yscale ("log")
    # ax2.set_ylim (0.01, 500)
    ax2.set_xlim (0, 1)
    ax2.plot (arm1, T1, lw=1, markersize=0, c="firebrick")
    ax2.plot (arm2, T2, lw=1, markersize=0, c="gold")
    fig2.savefig ("spinplots",dpi=1200, bbox_inches="tight")
    print (f"Plotted!")
    stop = time.time()
    print (f"Elapsed time is {stop-start} seconds.")



