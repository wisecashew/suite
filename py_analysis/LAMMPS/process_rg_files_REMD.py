#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import re
import os
import argparse

parser = argparse.ArgumentParser(description="Make the net radius of gyration files.")
args = parser.parse_args ()

if __name__=="__main__":

	TIMESTEP = -1
	DELTA    = 5000
	
	files = os.listdir(".")
	rg_files = list()
	for file in files:
		if re.match(r"^rad_gyr\.[0-9]+\.avg", file):
			rg_files.append(file)
	rg_files.sort(key=lambda x:int(x.split(".")[1]))

	print(f"Rg files to process: {rg_files}")

	for idx, rg_file in enumerate(rg_files):
		f = open(rg_file, 'r')
		new_file = open(f"TOT_RG_FILES/rg.{idx}.dump", 'w')
		hash_counter = 0
		for line in f:
			if re.match("^# Time-averaged", line):
				continue
			elif re.match("^# TimeStep", line):
				hash_counter += 1
				continue
			else:
				LINE = line.strip().split()
				LINE = [int(LINE[0]), float(LINE[1]), float(LINE[2]), float(LINE[3])]
				if hash_counter > 1:
					if LINE[0] == 0:
						continue
					else:
						LINE[0] = TIMESTEP + DELTA
				new_file.write(f"{LINE[0]} {LINE[1]} {LINE[2]} {LINE[3]}\n")
				TIMESTEP = LINE[0]
		new_file.close()
		f.close()
	
	

