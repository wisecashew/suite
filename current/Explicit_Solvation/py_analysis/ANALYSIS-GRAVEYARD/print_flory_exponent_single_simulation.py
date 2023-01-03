#!/home/satyend/.conda/envs/data_analysis/bin/python

import argparse 
import aux 
import scipy 
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt 
from sklearn.linear_model import LinearRegression 

parser = argparse.ArgumentParser (description="Print out the ensembled average of intermonomer distances, over a range of indices.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store',          help='enter a degree of polymerization.')
parser.add_argument('-U', metavar='U', dest='U', type=str, action='store',                help='enter potential energy surface.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store',                help='start parsing after this move number (not index or line number in file).')
parser.add_argument('-num', metavar='num', type=int, dest='num', action='store',          help='choose a trajectory.')
parser.add_argument('-T', metavar='T', dest='T', type=float, action='store',              help='choose a temperature.')
parser.add_argument('--coords', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('-d1', dest='d1', metavar='d1', action='store', type=int, help='starting difference index.') 
parser.add_argument('-d2', dest='d2', metavar='d2', action='store', type=int, help='ending difference index.'  )

args = parser.parse_args() 

def func (x, a, nu):
    return a * (x)**(nu)

x   = list(np.arange(args.d1, args.d2))
y   = [] 
U   = args.U
T   = args.T
num = args.num
dop = args.dop
coords_file    = args.c
starting_index = args.s
for i in x:
    y.append( aux.single_sim_flory_exp ( U, T, num, dop, coords_file, starting_index, i ) )

x = np.asarray (np.log(x)).reshape ((-1,1))
y = np.asarray (np.log(y))

print (y)

model = LinearRegression() 
model.fit(x, y)
r2 = model.score (x, y)
print (f"intercept: {model.intercept_}") 
print (f"slope:     {model.coef_}") 
print (f"R2:        {r2}")
