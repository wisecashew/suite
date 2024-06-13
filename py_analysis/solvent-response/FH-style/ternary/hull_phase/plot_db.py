import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpltern
import ternary
from scipy.optimize import fsolve
import FreeEnergy as FE
import pickle
import argparse 

parser = argparse.ArgumentParser(description='')
parser.add_argument('--db',  metavar='db',  dest='db',  type=str, action='store', help='location of database.')
parser.add_argument('--img', metavar='img', dest='img', type=str, action='store', help='name of image to be created (default: my_database).', default="my_database")
args = parser.parse_args()


if __name__=="__main__":
	df = pd.read_csv(args.db, sep='\|', engine='python', names=["vs", "vc", "vp", "chisc", "chips", "chipc", "phi_s", "phi_p", "l0", "l1", "l2", "phi_s1", "phi_p1", "phi_s2", "phi_p2", "phi_s3", "phi_p3", "w1", "w2", "w3"], skiprows=1)

	fig = plt.figure(figsize=(10,10))
	ax  = fig.add_subplot(projection="ternary")

	ps  = np.array(df["phi_s"].values,  dtype=float)
	pp  = np.array(df["phi_p"].values,  dtype=float)
	ps1 = np.array(df["phi_s1"].values, dtype=float)
	pp1 = np.array(df["phi_p1"].values, dtype=float)
	ps2 = np.array(df["phi_s2"].values, dtype=float)
	pp2 = np.array(df["phi_p2"].values, dtype=float)
	ps3 = np.array(df["phi_s3"].values, dtype=float)
	pp3 = np.array(df["phi_p3"].values, dtype=float)

	L = len(df["vs"].values)

	ax.scatter([], [], [], s=0.5, c='darkblue',   label="Single phase")
	ax.scatter([], [], [], s=0.5, c='coral',      label="Double phase")
	ax.scatter([], [], [], s=0.5, c='lightgreen', label="Triple phase")

	for i in range(L):
		if i % 100 == 0:
			print(f"Plotted  {i}th point out {L} points...", flush=True)
		if int(df["l0"].values[i]) == 1:
			ax.scatter(ps[i], 1-ps[i]-pp[i], pp[i], c='darkblue', s=0.8)
		elif int(df["l1"].values[i]) == 1:
			ax.scatter(ps[i], 1-ps[i]-pp[i], pp[i], c='coral',  s=0.8)
		elif int(df["l2"].values[i]) == 1:
			ax.scatter(ps[i], 1-ps[i]-pp[i], pp[i], c='lightgreen', s=0.8)
	ax.legend(loc='upper right', )
	fig.savefig(args.img, dpi=1000)
