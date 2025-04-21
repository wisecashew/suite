#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import pandas as pd
import os
import re
import argparse

parser = argparse.ArgumentParser(description="Make the net output files.")
args = parser.parse_args ()

if __name__=="__main__":

	TIMESTEP = -1
	DELTA    = 1000
	
	files        = os.listdir(".")
	output_files = list()
	for file in files:
		if re.match(r"^output.[0-9]+", file):
			output_files.append(file)
	output_files.sort(key=lambda x:int(x.split(".")[1]))
	
	print(f"Output files to process: {output_files}.")

	new_file = open(f"TOT_OUTPUT_FILES/output.dump", 'w')
	for idx,output_file in enumerate(output_files):
		f = open(output_file, 'r')
		for line in f:
			if re.match(r"^[0-9]+", line):
				LINE = [int(i) for i in line.strip().split()]
				if idx > 0:
					if LINE[0] == 0:
						continue
					elif idx > 0:
						LINE[0] = TIMESTEP + DELTA
				new_file.write(f"{LINE[0]}")
				for element in LINE[1:]:
					new_file.write(f" {element}")
				new_file.write("\n")
				TIMESTEP = LINE[0]
		f.close()
	new_file.close()
	
