#!/home/satyend/.conda/envs/data_analysis/bin/python

import pandas as pd 
import numpy as np 
import re
import aux
import argparse 

parser = argparse.ArgumentParser (description="Dump magnetization.")
parser.add_argument ("--file-path", dest='fp', action='store', type=str, help="Enter a file path to orientation file.")

args = parser.parse_args()

if __name__=="__main__":

	aux.dump_magnetization (args.fp, 1, 0)

