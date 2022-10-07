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
import itertools
from sklearn.linear_model import LinearRegression 

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

sys.stdout.flush() 

''' 
shebang for cluster: #!/usr/licensed/anaconda3/2020.7/bin/python
shebang for homemachine: #!/usr/bin/env python3
'''


import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the flory exponent from that file.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-U', metavar='U', dest='U', type=str, action='store', help='enter a degree of polymerization.')
parser.add_argument('-T', metavar='T', dest='T', type=float, action='store', help='enter a degree of polymerization.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this move number (not index or line number in file).', default=100)
parser.add_argument('--coords', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('-d1', dest='d1', metavar='d1', action='store', type=int, help='Starting index.')
parser.add_argument('-d2', dest='d2', metavar='d2', action='store', type=int, help='End index.')
parser.add_argument('--dump-file', dest='df', metavar='df', action='store', type=str, help='Name of dump file.')
args = parser.parse_args() 

divnorm = matplotlib.colors.SymLogNorm (0.005, vmin=-0.1, vmax=0.1)


def get_flory (U, T, num, dop, coords_file, starting_index, d1, d2):
    x = list (np.arange(d1, d2))
    y = [] 
    for i in x:
        y.append ( aux.single_sim_flory_exp (U, T, num, dop, coords_file, starting_index, i ) )

    x = np.asarray (np.log(x)).reshape((-1,1)) 
    y = np.asarray (np.log(y)) 

    model = LinearRegression() 
    model.fit(x, y)
    r2 = model.score (x, y) 
    return (model.coef_, r2) 


if __name__ == "__main__":    

    start = time.time() 
    ##################################

    U_list = [args.U]# aux.dir2U ( os.listdir (".") )
    PLOT_DICT = {} 
    dop            = args.dop
    coords_files   = args.c
    starting_index = args.s
    
    ######
    fig = plt.figure( figsize=(8,6) )
    ax  = plt.axes() 
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    i = 0 
    Tmax = [] 

    # instantiating pool
    pool1 = multiprocessing.Pool ( processes=50 )
    pool2 = multiprocessing.Pool ( processes=5  )
    
    pool_list = [pool1, pool2]
    
    f = open(args.df, "a") 

    for U in U_list:
        f.write ( "U = " + str(U) + ":\n" )
        print("Inside U = " + U + ", and N = " + str(dop), flush=True )
        # temperatures = aux.dir2float ( os.listdir( str(U) +"/DOP_"+str(dop) ) )
        temperatures = [args.T]
        # get num_list for each temperature 
        master_temp_list = []
        master_num_list = []
        flory_dict    = {}
        r2_dict       = {} 
        ntraj_dict = {}
        flory_mean = [] 
        flory_r2   = [] 
        for T in temperatures: 
            # print ("T is " + str(T), flush=True) 
            num_list = list(np.unique ( aux.dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
            master_num_list.extend ( num_list )
            master_temp_list.extend ( [T]*len( num_list ) )
            ntraj_dict[T] = len ( num_list )
            flory_dict[T] = []
            r2_dict   [T] = [] 


        # start multiprocessing... keeping in mind that each node only has 96 cores 
        # start splitting up master_num_list and master_temp_list 
        mtemp_list_p1 = master_temp_list[0:50] 
        mtemp_list_p2 = master_temp_list[50:100]
        mtemp_list_p3 = master_temp_list[100:105]
        mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3] 

        mnum_list_p1  = master_num_list[0:50]
        mnum_list_p2  = master_num_list[50:100]
        mnum_list_p3  = master_num_list[100:105] 
        mnum_list     = [mnum_list_p1, mnum_list_p2, mnum_list_p3] 

        # it is a shitty dict 
        shitty_dict = {0:0, 1:0, 2:1 }

        for uidx in range(3):
            
            # pool = multiprocessing.Pool ( processes=len(mtemp_list[uidx]) ) # len(num_list)) 
            #with pool as p:
            results = pool_list[ shitty_dict[uidx] ] .starmap ( get_flory, zip( itertools.repeat(U), mtemp_list[uidx],\
                     mnum_list[uidx], itertools.repeat(dop), \
                    itertools.repeat(coords_files), itertools.repeat(starting_index), itertools.repeat(args.d1), itertools.repeat(args.d2) ) )

            # print (len(results))
            # pool.join() 

            print ("Pool has been closed. This pool has {} threads.".format (len(results) ), flush=True )     

            for k in range(len(mtemp_list[uidx])):
                flory_dict[mtemp_list[uidx][k]].append( results[k][0] )
                r2_dict[mtemp_list[uidx][k]].append   ( results[k][1] ) 
        
            for T in np.unique (mtemp_list[uidx]):
                flory_mean.append( np.mean ( flory_dict[T] ) ) 
                flory_r2.append (np.mean( r2_dict[T] ) )

        PLOT_DICT [U] = (np.asarray(flory_mean), np.asarray(flory_r2))
        
        f.write ("flory exp: ")
        for elem in flory_mean:
            f.write ( "{: 2.2f} ".format(elem))
        f.write ("\n")
        f.write ("r2:        ")
        for elem in flory_r2:
            f.write ( "{: 2.2f} ".format(elem) )
        f.write ("\n" )
        f.write ("T:         ")
        for elem in temperatures: 
            f.write ( "{: 2.2f} ".format(elem) ) 
        f.write("\n") 
        i+=1 
        f.flush()

    pool1.close()
    pool1.join()

    pool2.close()
    pool2.join()

    f.close()
    
    i=0



    ##################################
    stop = time.time() 

    print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)
        
