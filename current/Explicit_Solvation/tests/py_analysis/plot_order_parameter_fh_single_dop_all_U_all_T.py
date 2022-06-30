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
# parser.add_argument('--excl-vol' , dest='ev', action ='store_true' , help='Flag to include excluded volume forcefield.', default=False) 
parser.add_argument('--ortn-file', dest='of', metavar='orientation', action='store', type=str, help='Name of orientation dump file to parse information.', default='orientation')
parser.add_argument('--show-plot', dest='sp', action ='store_true' , help='Flag to include to see plot.') 
args = parser.parse_args() 


if __name__=="__main__":

	start = time.time() 
	#######################################

	fig = plt.figure( figsize=(8,6) )
	ax  = plt.axes  () 

	U_list    = aux.dir2U ( os.listdir (".") ) 
	PLOT_DICT = {} 

	starting_index = args.s
	ortn_file      = args.of
	i    = 0
	Tmax = 0 

	# instantiating pool 
	pool1 = multiprocessing.Pool ( processes=50 ) 
	pool2 = multiprocessing.Pool ( processes=5  )

	pool_list = [ pool1, pool2 ]

	for U in U_list:

		temperatures = aux.dir2float ( os.listdir ( str(U) + "/DOP_" + str(dop) ) )
		# get the num_list for each temperature 
		master_temp_list = [] 
		master_num_list  = [] 
		ord_par_dict     = {} 
		ntraj_dict       = {} 

		# define stores 
		ord_par_mean     = [] 
		ord_par_std      = [] 

		for T in temperatures:
			num_list         = list(np.unique ( dir2nsim (os.listdir (str(U) + "/DOP_" + str(dop) + "/" + str(T) ) ) ) )
			master_num_list.extend  ( num_list )
			master_temp_list.extend ( [T]* len(num_list) )
			ntraj_dict[T]    = len (num_list)
			ord_par_dict [T] = []

		# start multiprocessing... keeping in mind that each node only has 96 cores
		# start splitting up master_num_list and master_temp_list 
		mtemp_list_p1 = master_temp_list [0:50]
		mtemp_list_p2 = master_temp_list [50:100]
		mtemp_list_p3 = master_temp_list [100:105]
		mtemp_list    = [mtemp_list_p1, mtemp_list_p2, mtemp_list_p3]

		mnum_list_p1  = master_num_list [0:50]
		mnum_list_p2  = master_num_list [50:100]
		mnum_list_p3  = master_num_list [100:105]

		# define a shitty dict 
		shitty_dict = {0:0, 1:0, 2:1}

		for uidx in range(3):

			results = pool_list [ shitty_dict[uidx] ].starmap ( aux.obtain_order_parameter, zip( itertools.repeat(U), itertools.repeat(dop), \
				mtemp_list[uidx], itertools.repeat(ortn_file), mnum_list[uidx], itertools.repeat(ortn_file) ) ) 

			print ("Pool has been closed. This pool has {} threads".format ( len(results ) ), flush=True )

			for k in range ( len (mtemp_list[uidx]) ):
				ord_par_dict[ mtemp_list[uidx][k] ].append ( results[k] )

			for T in np.unique ( mtemp_list [uidx] ):
				ord_par_mean.append ( np.mean ( ord_par_dict[T] ) ) 
				ord_par_std.append  ( np.std  ( ord_par_dict[T] ) / np.sqrt ( ntraj_dict[T] ) )

		PLOT_DICT [U] = (np.asarray (ord_par_mean), np.asarray (ord_par_std))

	pool1.close()
	pool1.join () 

	pool2.close() 
	pool2.join ()

	i = 0 

	for U in U_list: 
		ax.errorbar ( temperatures, PLOT_DICT[U][0], yerr=PLOT_DICT[U][1], fmt='o', markeredgecolor='k', \
			linestyle='-', elinewidth=1, capsize=0, linewidth=1, color=cm.seismic(i/len(U_list)), label='_nolegend_' )
		i += 1

	##############################################################

	my_cmap = cm.seismic
	sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )
    cbar = plt.colorbar(sm, orientation='vertical') 
    cbar.set_ticks ( [0, 1] )
    cbar.set_ticklabels( ["Weakest", "Strongest"] ) 
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel ("Strength of aligned \nmonomer-solvent interactions", fontsize=18, rotation=270)
    ax.set_xscale('log')
    ax.set_xlabel ( "Temperature (reduced)", fontsize=18) 
    ax.set_ylabel ( "$\\langle R_g^2 \\rangle/ \\langle R_g ^2 \\rangle _{\\mathrm{max}}$", fontsize=18)     
    ax.set_yticks (np.linspace(0, 1, 11)) 
    plt.savefig   ( "DOP_"+str(dop)+"_multiple_rg.png", dpi=1000)
    
    if show_plot_bool:
        plt.show() 

    ##################################
    stop = time.time() 
	
	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)



