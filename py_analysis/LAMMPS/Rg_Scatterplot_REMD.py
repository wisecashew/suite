#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the Rg histogram.")
parser.add_argument('--image',     dest='img',    action='store', type=str,   help='Name of image to be created.', default=None)
args = parser.parse_args()


if __name__=="__main__":

	assert args.img != None

	# track how the temperature index goes around
	# indices go from 0 to 39
	T_range   = np.linspace(275, 360, 40)
	RG_over_T = dict()
	for idx in range(0, 40):
		RG_over_T[idx] = [] # list()

	DF_list = [] # list()
	for idx in range(0, 40):
		df = pd.read_csv("TOT_RG_FILES/rg."+str(idx)+".dat", sep='\s+', names=["timestep", "rg", "rgx", "rgy", "rgz"], engine="python", skiprows=2000)
		DF_list.append(df)
	
	print(f"Temperature = ", end='', flush=True)
	for Tidx in range(0, 40):
		print(f"{T_range[Tidx]:.2f}, ",end='', flush=True)
		output_file = open("TOT_OUTPUT_FILES/output.net", 'r')
		for line in output_file:
			if re.match(r"^[0-9]+", line):
				line       = line.strip().split()
				line       = [int(i) for i in line]
				ts         = line[0]
				indices    = line[1:]
				replica_id = indices.index(Tidx)
				timesteps  = np.array(DF_list[replica_id]["timestep"].values, dtype=int)
				if np.sum(timesteps==ts) == 1:
					RG_over_T[Tidx].extend(DF_list[replica_id]["rg"].values[timesteps==ts].tolist())
		output_file.close()

	print()
	print(f"Make images...", flush=True)

	for Tidx in range(0, 40):
		print(f"@ T = {T_range[Tidx]:.2f} and ", end='')
		fig = plt.figure(num=Tidx, figsize=(3,3))
		ax  = plt.axes()
		print(f"number of samples = {len(RG_over_T[Tidx])}.")
		ax.scatter(range(len(RG_over_T[Tidx])), RG_over_T[Tidx], c='steelblue', marker='o', alpha=0.5)
		ax.axhline(y=12.5, ls="--", c="darkred")
		ax.set_title(f"T = {T_range[Tidx]:.2f} K")
		ax.set_xlabel("Timesteps (every 5ps)")
		ax.set_ylabel("Rg")
		ax.set_ylim(8, 20)
		fig.savefig("TOT_IMAGES/"+args.img+f"_{T_range[Tidx]:.2f}.png", dpi=1200, bbox_inches="tight")
		plt.close(fig)
