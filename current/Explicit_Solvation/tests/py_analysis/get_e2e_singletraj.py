#!/home/satyend/.conda/envs/data_analysis/bin/python

import time 
import sys 
import aux 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''


import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-U', metavar='UX', dest='U', type=str, action='store', help='Enter a potential energy surface.')
parser.add_argument('-T', metavar='T', dest='T', type=float, action='store', help='Enter a temperature.') 
parser.add_argument('-n', metavar='N', dest='n', type=int, action='store', help='Enter a trajectory number.') 
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this index.', default=100)
parser.add_argument('--coords', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('--show-plot', dest='sp', action='store_true', help='Flag to include to see plot.') 
args = parser.parse_args() 


if __name__ == "__main__":    

    start = time.time() 

    aux.get_e2e_data ( args.U, args.T, args.dop, args.n, args.s, args.c, args.sp )

    stop = time.time() 

    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)
        
