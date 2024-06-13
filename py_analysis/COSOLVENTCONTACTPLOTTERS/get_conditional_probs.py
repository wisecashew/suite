#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd

def count_occurrences(condensed):
	counts    = dict()
	for entry in condensed:
		x1, x2, x3 = entry[0], entry[1], entry[2]
		counts[(x1, x2, x3)] = counts.get((x1, x2, x3), 0) + 1
	return counts

def calculate_joint_probs(counts, total_entries):
	joint_probs = {}
	for (x1, x2, x3), count in counts.items():
		joint_probs[(x1, x2)] = joint_probs.get((x1, x2), 0) + count / total_entries
	return joint_probs

def calculate_cond_probs(counts, joint_probs):
	cond_probs = {}
	for (x1, x2, x3), count in counts.items():
		cond_probs[(x1, x2, x3)] = count/joint_probs[(x1, x2)]
	return cond_probs

if __name__=="__main__":

	df            = pd.read_csv("SOLVATION.csv", sep='\|', names=["U", "bulkfrac", "N_s", "T_ms", "N_ms"], skiprows=1, engine='python')
	condensed     = np.array([df["N_s"].values, df["T_ms"].values, df["N_ms"].values], dtype=int).T
	total_entries = len(condensed)

	counts            = count_occurrences(condensed)
	joint_probs       = calculate_joint_probs(counts, total_entries)
	conditional_probs = calculate_cond_probs (counts, joint_probs)

	RN_cond_probs = dict()
	normalization = dict()
	Nms_range = np.linspace(0, 771, 1, dtype=int)
	for key in joint_probs:
		normalization[key] = 0
		for Nms in Nms_range:
			RN_key = (key[0], key[1], Nms)
			RN_cond_probs[RN_key] = conditional_probs.get(RN_key, 0)
			normalization[key]   += RN_cond_probs[RN_key] 
	
	for key in joint_probs:
		normalization[key] = 0
		for Nms in Nms_range:
			RN_key = (key[0], key[1], Nms)
			RN_cond_probs[RN_key] = RN_cond_probs[RN_key]/normalization[key]


