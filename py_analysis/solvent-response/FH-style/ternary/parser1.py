import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpltern

if __name__=="__main__":

	lsize = 3
	plt.rcParams['font.family'] = 'Arial'
	font = {'color':  'black',
        'weight': 'normal',
        'size': lsize}

	fig = plt.figure(num=1, figsize=(5,5))
	ax  = fig.add_subplot (projection="ternary")
	
	df = pd.read_csv ("comparing_mu.txt", sep='\s+', engine="python", skiprows=1, names=["i1", "i2", "dmu", "mu_a1", "mu_b1", "mu_c1", \
		"phi_a1", "phi_b1", "phi_c1", "mu_a2", "mu_b2", "mu_c2", "phi_a2", "phi_b2", "phi_c2"])

	print (df)
	df = df.loc[df["dmu"]<0.01]
	print (df)

	phi_a = df["phi_a1"].values; phi_an = df["phi_a2"].values
	phi_b = df["phi_b1"].values; phi_bn = df["phi_b2"].values
	phi_c = df["phi_c1"].values; phi_cn = df["phi_c2"].values

	ax.scatter (phi_a, phi_b, phi_c, s=0.5, c="steelblue")

	for i in range (len (phi_a)):
		p1 = [phi_a[i], phi_an[i]] 
		p2 = [phi_b[i], phi_bn[i]]
		p3 = [phi_c[i], phi_cn[i]]
		ax.plot (p1, p2, p3, lw=1, ls="--", c='k')


	ax.set_tlabel('Vol. frac. A')
	ax.set_llabel('Vol. frac. B')
	ax.set_rlabel('Vol. frac. C')

	# Set axis limits
	ax.set_tlim(0, 1)
	ax.set_llim(0, 1)
	ax.set_rlim(0, 1)
	# ax.ticks(axis='lbr', multiple=5, linewidth=1, offset=0.025)

	positions = ['tick1', 'tick2']
	for position in positions:
	    ax.taxis.set_ticks_position(position)
	    ax.laxis.set_ticks_position(position)
	    ax.raxis.set_ticks_position(position)

	ax.grid()

	plt.savefig ("estimated_binodal", dpi=1200)
