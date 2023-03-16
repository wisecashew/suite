#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np
import pandas as pd 
import shutil 
import argparse 
import aux 
import os 


parser = argparse.ArgumentParser(description="Go into the low temperature simulations and only keep the copies of the least energetic simulation.")
parser.add_argument('-dop', dest='dop', action='store', type=int, help='Provide degree of polymerization.') 
parser.add_argument('-U', dest='U', action='store', type=str, help='Provide potential energy surface.')
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')

args = parser.parse_args() 

if __name__=="__main__": 

	# get the entire list of potential energy surfaces 
	U_list = [args.U] 

	# instantiate plt figure 
	# plt.figure ( figsize=(8,6) )
	temperatures = [0.01, 0.05, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
	name_list    = ["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "ms1s2_tot",  "ms1s2_aligned", "ms1s2_naligned", "time_step"]
	for U in U_list:
		min_energy = +1
		for temp in temperatures:
			# find the trajectory with the most contacts 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(args.dop) + "/" + str(temp) ) ) )

			for num in num_list: 
				df = pd.read_csv( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc", sep=' \| ', names = name_list, engine='python', skiprows=args.s) 
				net_energy = df["energy"].values[-1]
				if  net_energy < min_energy:
					min_energy   = net_energy
					energy_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc"
					coords_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/coords_"+str(num)+".mc"
					lattice_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/lattice_dump_"+str(num)+".mc"
			new_energy = 0
			for num in num_list:

				destination_energy_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)+".mc"
				destination_coords_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/coords_"+str(num)+".mc"
				destination_lattice_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/lattice_dump_"+str(num)+".mc" 

				if energy_dump == destination_energy_dump or coords_dump == destination_coords_dump or lattice_dump == destination_lattice_dump:
					continue
				else:
					try:
						shutil.copyfile(energy_dump,  destination_energy_dump)
						shutil.copyfile(coords_dump,  destination_coords_dump)
						shutil.copyfile(lattice_dump, destination_lattice_dump)
						print("File copied successfully.")
				 
					# If source and destination are same
					except shutil.SameFileError:
				    	    print("Source and destination represents the same file.")
				 
					# If destination is a directory.
					except IsADirectoryError:
					    print("Destination is a directory.")
					 
					# If there is any permission issue
					except PermissionError:
					    print("Permission denied.")
					 
					# For other errors
					except:
					    print("Error occurred while copying file.")
