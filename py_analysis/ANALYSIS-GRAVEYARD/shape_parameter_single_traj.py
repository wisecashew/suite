#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np 
import re 
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib.cm as cm 
import os
import argparse 
import aux
import multiprocessing
import itertools
import sys

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the average shape parameter for this trajectory.")
parser.add_argument('-U', metavar='U', dest='U', type=str, action='store', help='Enter a potential energy surface.' )
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-T', metavar='T', dest='T', type=float, action='store', help='Enter a temperature')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number.', default=100)
parser.add_argument('--coords-file', dest='e', action='store', type=str, metavar='coords', help='Name of coords file to parse through.') 
parser.add_argument('-num', dest='num', action='store', type=int, metavar='N', help='Trajectory number.')
# parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include if you want the image to be rendered on a screen.', default=False)

args = parser.parse_args() 


def get_shape_parameter_from_one_traj ( U, DOP, T, num, coords_file, starting_index ):
    
    xlen = aux.edge_length (DOP) 
    ylen = aux.edge_length (DOP)
    zlen = aux.edge_length (DOP)

    filename = U + "/DOP_" + str(DOP) + "/" + str(T) + "/" + coords_file  + "_" + str(num) 
    master_dict = aux.get_pdict (filename, starting_index, DOP, xlen, ylen, zlen)
    
    # for key in master_dict:
    #     print(key)

    shape_parameter = aux.get_shape_param ( master_dict, xlen, ylen, zlen ) 

    return shape_parameter 

##############################################################################################


if __name__=="__main__":
    
    shape_param = get_shape_parameter_from_one_traj ( args.U, args.dop, args.T, args.num, args.e, args.s ) 

    print ("Shape parameter for U = " + args.U +  ",  T = " + str(args.T) + ", ntraj = " + str(args.num) + " is ", shape_param) 

