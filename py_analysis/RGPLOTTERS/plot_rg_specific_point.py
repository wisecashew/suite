#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np 
import re 
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import os
import sys 
sys.path.insert(0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux 
import time 
import multiprocessing 
import itertools

# os.system("taskset -p 0xfffff %d" % os.getpid())
# os.environ['MKL_NUM_THREADS'] = '1'
# os.environ['NUMEXPR_NUM_THREADS'] = '1'
# os.environ['OMP_NUM_THREADS'] = '1'
# sys.stdout.flush()

'''
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''

import argparse
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-s',             metavar='S',    type=int,              dest='s',       action='store', help='Calculate Rg at this point.', default=100)
parser.add_argument('--dop',          metavar='dop',  type=int,              dest='dop',     action='store', help='Degree of polymerization.', default=None)
parser.add_argument('--coords',       dest='c',       metavar='coords.txt',  action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
args = parser.parse_args() 

def infiltrate_coords_get_rg (coords_file, index, dop):

	filename    = coords_file
	edge        = aux.edge_length (dop)
	master_dict = aux.get_pdict (filename, index, dop, edge, edge, edge)
	coord_arr   = aux.unfuck_polymer(master_dict[index][0], edge, edge, edge)
	# print(coord_arr)
	r_com       = np.mean(coord_arr, axis=0) # get center of mass 
	offset      = coord_arr - r_com
	rg          = np.sqrt(np.sum (np.square (offset) )/ dop) # added the np.sqrt
	return rg 



if __name__ == "__main__":

##################################
	dop          = args.dop
	coords_files = args.c
	index        = args.s

	rg = infiltrate_coords_get_rg(coords_files, index, dop)
	print(f"Radius of gyration is {rg}.", flush=True)
######
