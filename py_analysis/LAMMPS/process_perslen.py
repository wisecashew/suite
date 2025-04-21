#!/home/satyend/.conda/envs/phase/bin/python

import re
import numpy as np
import pandas as pd
import argparse
import pickle
import time
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser(description="Get the information necessary for persistence length.")
parser.add_argument("--pl-dir", dest="pld", type=str, action="store", required=True, help="enter address of directory with pkl files. make sure it doesn't end with '/'.")
args = parser.parse_args()

def fit_exponential(x_data, y_data, lb):
	def model(x, lp):
		return np.exp(-x*lb/lp)

	# --- set up some initial guesses ---
	lp_init = np.mean(x_data)	# mean of x-data is a guess

	# --- perform the curve fitting --- 
	lp_opt, cov = curve_fit(model, x_data, y_data, p0=[lp_init])

	return lp_opt[0], cov

if __name__=="__main__":
	
	# ~~~ set up the timer ~~~
	start     = time.time() 
	max_delta = 39

	# --- get the list of files ---
	print(f"Get the list of files.", flush=True)
	list_of_files = os.listdir(args.pld)	# get the list of files
	pl_databases   = list()					# a place to store all the relevant files
	for file in list_of_files:				# loop through all the files in the directory
		if file.endswith(".pkl"):			# check if a file ends with "pkl"
			pl_databases.append(file)		# adding all pkl files to pl_databases

	# --- setting up the structure for the master pl databse ---
	print(f"Set up the structure of the ultimate database.", flush=True)
	T_range = np.linspace(275, 360, 40)			# get all the relevant temperatures
	master_pl_db = dict()						# instantiate the object
	for Tidx in range(40):						# loop through each temperature index
		master_pl_db[Tidx] = dict()				# seting up the structure...
		master_pl_db[Tidx]["bond_length"] = []	# instantiate the bond_length container
		for delta in range(1,max_delta):				# loop through each delta
			master_pl_db[Tidx][delta] = []		# set up the structure

	# --- adding data to the data structure ---
	print(f"Populate the ultimate database.", flush=True)
	for pl_db in pl_databases:						# loop through all the database files
		my_db    = open(args.pld+"/"+pl_db, 'rb')	# read the database file
		local_db = pickle.load(my_db)				# load up the database
		my_db.close()								# close the database file
		for Tidx in range(40):						# loop through each temperature index
			master_pl_db[Tidx]["bond_length"].extend(local_db[Tidx]["bond_length"])		# add the average bond lengths at each temp to a master ds
			for delta in range(1, max_delta):													# loop through each delta
				master_pl_db[Tidx][delta].extend(local_db[Tidx][delta])					# set up the database

	# --- process the data structure to get correlations
	print(f"Appropriately dissect the database and make plots.", flush=True)
	ave_bond_length = dict()	# instantiate a dictionary to store average bond lengths at each temperature
	correlations    = dict()	# instatiate a dictionary to store correlations at each temperature
	for Tidx in range(40):		# loop through all temperature indices
		ave_bond_length[Tidx] = np.mean(master_pl_db[Tidx]["bond_length"])		# add the average bond length to the dictionary
		correlations[Tidx]    = []												# instantiate the correlations dictionary
		for delta in range(0, max_delta):												# loop through all the deltas
			if delta == 0:														# delta = 0 is a special case (correlation = 1)
				correlations[Tidx].append(1)									# setting correlation to 1 for delta = 0
			else:																# for delta > 0
				correlations[Tidx].append(np.mean(master_pl_db[Tidx][delta]))	# set the correlation to the mean of the correlations calculated
		lb     = ave_bond_length[Tidx]
		deltas = np.array(range(max_delta))
		corrs  = np.array(correlations[Tidx])
		lp, cov = fit_exponential(deltas, corrs, ave_bond_length[Tidx])
		print(f"The average bond length for T = {T_range[Tidx]:.2f} K is {ave_bond_length[Tidx]:.2f} and... ", end='', flush=True)
		print(f"the persistence length is {lp}.", flush=True)

		# --- plot the correlations ---
		print(f"Making image for T = {T_range[Tidx]:.2f} K...", end=' ', flush=True)
		fig = plt.figure(figsize=(3,3), num=Tidx)	# set up the figure with a size and number
		ax  = plt.axes()							# set up the axes
		ax.plot(deltas, corrs, color='steelblue', \
		mec='k', marker='o', lw=1, markersize=4, ls='--', label=f"Data")								# plot the correlations
		ax.plot(deltas, np.exp(-deltas*lb/lp), lw=2, color='darkred', alpha=0.5, label=f"Fit $l_p$ = {lp:.2f}, $l_b$ = {lb:.2f}")					# plot the exponential fit
		ax.legend(loc="upper right", prop={"size": 6})
		ax.set_xlabel("$\\delta$")
		ax.set_ylabel("C($\\delta$)")
		ax.set_title(f"T = {T_range[Tidx]:.2f}")
		ax.set_xlim(0, 40)
		ax.set_ylim(np.min(corrs), 1)
		ax.set_xticks([0, 10, 20, 30, 40])
		if np.min(corrs)<0:
			ax.set_yticks([np.min(corrs), 0.0, 0.2 ,0.4, 0.6, 0.8, 1.0])
		else:
			ax.set_yticks([0.0, 0.2 ,0.4, 0.6, 0.8, 1.0])

		fig.savefig(args.pld+"/"+f"corr_{Tidx}.png", dpi=1200, bbox_inches="tight")		# save the figure
		plt.close(fig)																	# close the figure
		print("done!", flush=True)

	# ~~~ get the end time ~~~
	stop = time.time()
	print(f"Completed computation. Time for computation is {stop-start} seconds.", flush=True)
