#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import argparse

if __name__=="__main__":

	f = open("TOT_HBONDS/hbonds.0.pkl", 'rb')
	data = pickle.load(f)
	f.close()

	hbond_data = data[0].results.hbonds

	for frame in np.unique(hbond_data[:,0]):
		frame_hbonds = hbond_data[hbond_data[:,0]==frame]
		print(f"frame = {frame}.")
		print(f"total = {frame_hbonds.shape[0]}")
		print(f"all info = {frame_hbonds}")
		exit()



