import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpltern
import ternary
from scipy.optimize import fsolve
import FreeEnergy as FE
import pickle

if __name__=="__main__":
	df = pd.read_csv("double.db", sep='\|', engine='python', names=["vs", "vc", "vp", "chisc", "chips", "chipc", "phi_s", "phi_p", "l0", "l1", "l2", "phi_s1", "phi_p1", "phi_s2", "phi_p2", "phi_s3", "phi_p3", "w1", "w2", "w3"], skiprows=1)

	fig = plt.figure(figsize=(3,3))
	ax  = fig.add_subplot()
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)

	p1 = np.array([df["phi_s1"].values, df["phi_p1"].values]).T
	p2 = np.array([df["phi_s2"].values, df["phi_p2"].values]).T
	p1, keep = ternary.remove_close_rows(p1)
	p2       = p2[keep]

	for 

	# for i in range(len(p1)):
	# 	ax.plot([p1[i,0], p2[i,0]], [p1[i,1], p2[i,1]], c='pink', ls='--', lw=0.5)
	ax.scatter(p1[:,0], p1[:,1], c="limegreen", s=0.1, alpha=0.1)
	ax.scatter(p2[:,0], p2[:,1], c="darkgreen",  s=0.1, alpha=0.5)
	# ax.scatter(p1[15:35,0], p1[15:35,1], c="gold", s=0.5)
	# ax.scatter(p2[15:35,0], p2[15:35,1], c="red",  s=0.5)


	fig.savefig("two_split.png", dpi=1000, bbox_inches="tight")
