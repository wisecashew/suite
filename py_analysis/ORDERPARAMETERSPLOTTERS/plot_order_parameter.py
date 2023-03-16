#!/home/satyend/.conda/envs/data_analysis/bin/python

import numpy as np 
import re 
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import pandas as pd
import os
import aux 
import time 
import sys 
import multiprocessing 
import cmath 

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-dop', metavar='DOP'   , dest='dop', type=int , action='store', help='enter a degree of polymerization.')
parser.add_argument('-s'  , metavar='S'     , type=int  , dest='s' , action='store', help='start parsing after this index.', default=100)
parser.add_argument('--excl-vol' , dest='ev', action ='store_true' , help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--ortn-file', dest='of', metavar='orientation', action='store', type=str, help='Name of orientation dump file to parse information.', default='orientation')
parser.add_argument('--show-plot', dest='sp', action ='store_true' , help='Flag to include to see plot.') 
args = parser.parse_args() 


if __name__=="__main__":

	pi = np.pi 

	# store the potential energy surface in a variable 
	U_list = aux.dir2U ( os.listdir(".") )
	N              = args.dop 
	ortn_file      = args.of
	starting_index = args.s

	start_str = "START for Step"
	end_str   = "END" 

	f = open("orientation_1", 'r')

	order_parameter_list = [] 
	extract_orr = False 
	start_bool  = False 
	for line in f:

		if re.match ( start_str, line ):
			a = re.search ("\d+", line)
			extract_orr = True
			if int( a.group(0) ) == starting_index:
				start_bool = True 
			oparam = 0 
			count  = 0 

		elif re.match ( end_str, line ) and start_bool:
			# print ( line )
			extract_orr = False 
			order_parameter_list.append( oparam/count ) 

		elif extract_orr and start_bool: 
			or_list = aux.extract_loc_from_string( line )

			for cnum in or_list:
				oparam += complex ( np.cos(2*pi*cnum/6 ), np.sin (2*pi*cnum/6) )
				count  += 1 

	f.close() 

	print ( "Order parameter for solvation shell is: ", abs(np.mean (order_parameter_list[400:]) ) ) 

			


			
