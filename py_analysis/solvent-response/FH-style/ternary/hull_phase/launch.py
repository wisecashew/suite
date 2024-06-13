import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import tangent
import mpltern
import pickle
import phase
import argparse

parser = argparse.ArgumentParser(description='Create a skeleton solution for the binodal. This is a memory-intensive computation.'         )
parser.add_argument('--chisc',               metavar='chi_sc',  dest='chi_sc',        type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',               metavar='chi_ps',  dest='chi_ps',        type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',               metavar='chi_pc',  dest='chi_pc',        type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',                   metavar='vs',      dest='vs',            type=float, action='store', help='specific volume of solvent.'  )
parser.add_argument('-vc',                   metavar='vc',      dest='vc',            type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',                   metavar='vp',      dest='vp',            type=float, action='store', help='specific volume of polymer.'  )
args = parser.parse_args() 

###################################################
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
###################################################

if __name__=="__main__":

    start = time.time()
    inputs = dict()
    inputs["chi_sc"] = 2.4 # args.chi_sc
    inputs["chi_ps"] = 2.4 # args.chi_ps
    inputs["chi_pc"] = 2.4 # args.chi_pc
    inputs["vs"]     = 1 # args.vs
    inputs["vp"]     = 1 # args.vp
    inputs["vc"]     = 1 # args.vc

    # set up objects and crit points
    print(f"Setting up the Phase object...", flush=True, end=' ')
    P = phase.Phase(inputs)
    print("done!", flush=True)

    df = pd.read_csv("double.db", sep='\|', engine='python', names=["vs", "vc", "vp", "chisc", "chips", "chipc", "phi_s", "phi_p", "l0", "l1", "l2", "phi_s1", "phi_p1", "phi_s2", "phi_p2", "phi_s3", "phi_p3", "w1", "w2", "w3"], skiprows=1)
    length = len(df["vs"].values)
    arm_1 = np.empty((0,2))
    arm_2 = np.empty((0,2))

    for i in range(length):
        def delta_mu(phi):
            eq1 = P.sym_mu_ps.delta_mu_s(df["phi_s1"].values[i], phi[0], phi[1], phi[2])
            eq2 = P.sym_mu_ps.delta_mu_p(df["phi_s1"].values[i], phi[0], phi[1], phi[2])
            eq3 = P.sym_mu_ps.delta_mu_c(df["phi_s1"].values[i], phi[0], phi[1], phi[2])
            return [eq1, eq2, eq3]
        root = fsolve(delta_mu, [df["phi_p1"].values[i], df["phi_s2"].values[i], df["phi_p2"].values[i]])
        p1   = np.array([df["phi_s1"].values[i], root[0]])
        p2   = np.array([root[1], root[2]])

        if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0 or root[2] > 1 or root[2] < 0:
            print(f"Excessive roots found for index {i}...", flush=True)
            continue
        elif (np.abs(delta_mu(root))>1e-6).any():
            print(f"Not good enough roots found for index {i}...", flush=True)
            continue
        elif np.linalg.norm(p1-p2) < 1e-6:
            print(f"Roots are too close...", flush=True)
            print(f"p1 = {p1}", flush=True)
            print(f"p2 = {p2}", flush=True)
        else:
            print(f"Found root!", flush=True)
            print(f"p1 = {p1}", flush=True)
            print(f"p2 = {p2}", flush=True)
            arm_1 = np.vstack((arm_1, p1))
            arm_2 = np.vstack((arm_2, p2))
    

    fig = plt.figure(figsize=(3,3))
    ax  = fig.add_subplot()
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    print(f"Length of points = {len(arm_1)}")
    ax.scatter(arm_1[:,0], arm_1[:,1], c="limegreen", s=0.1)
    ax.scatter(arm_2[:,0], arm_2[:,1], c="darkgreen",  s=0.1)
    for i in range(len(arm_1)):
        ax.plot([arm_1[i,0], arm_2[i,0]], [arm_1[i,1], arm_2[i,1]], lw=0.5, c='violet')

    fig.savefig("diagram.png", dpi=1000, bbox_inches="tight")





