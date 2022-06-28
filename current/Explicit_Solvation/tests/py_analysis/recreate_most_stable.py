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
parser.add_argument('-T', dest='T', action='store', type=float, help='Provide a temperature.')
parser.add_argument('-s', dest='s', action='store', type=int, help='Provide a starting index from when to sample.', default=100)
parser.add_argument('--dump-file', dest='e', metavar='energydump', action='store', type=str, help='Name of energy dump file to parse information.', default='energydump.txt')
parser.add_argument('--sim-type', dest='st', metavar='entropy', action='store', type=str, help='Enter the type of simulation that is being run.')

args = parser.parse_args() 

if args.st == "entropy":
	name_list = [ "energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step" ]

elif args.st == "FH": 
	name_list = [ "energy", "mm_tot", "ms_tot", "time_step" ]

elif args.st == "cosolvent":
	name_list = [ "energy", "mm_tot", "mm_aligned", "mm_naligned", "ms1_tot", "ms1_aligned", "ms1_naligned", "ms2_tot", "ms2_aligned", "ms2_naligned", "time_step" ]

elif args.st == "comonomer": 
	name_list = [ "energy", "mm_tot", "m1m1", "m2m2", "m1m2", "ms_tot", "m1s1", "m2s1", "time_step" ] 

else :
	print ("simulation type not found. Exiting...")
	exit() 


if __name__=="__main__": 

	# get the entire list of potential energy surfaces 
	U_list = [args.U] 

	# instantiate plt figure 
	# plt.figure ( figsize=(8,6) )

	for U in U_list:
		min_energy = +1
		temperatures = np.unique ( aux.dir2float( os.listdir (str(U) + "/DOP_" + str(args.dop) ) ) )
		temperatures = [args.T]
		for temp in temperatures:
			# find the trajectory with the most contacts 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(args.dop) + "/" + str(temp) ) ) )

			for num in num_list: 
				df = pd.read_csv( str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num), sep=' \| ', names = name_list, engine='python', skiprows=args.s) 
				net_energy = np.min( df["energy"].values ) 
				if  net_energy < min_energy:
					min_energy = net_energy
					energy_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)
					coords_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/coords_"+str(num)
					lattice_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/lattice_dump_"+str(num)


	print ("Minimum energy found in \"" + energy_dump +"\"")


	# get the coords, and copy them over 

	source_energy_dump   = energy_dump
	source_coords_dump   = coords_dump
	source_lattice_dump  = lattice_dump

	for U in U_list:
		for temp in temperatures: 
			num_list = np.unique ( aux.dir2nsim ( os.listdir ( str(U) + "/DOP_" + str(args.dop) + "/" + str(temp) ) ) )

			for num in num_list: 

				destination_energy_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/"+args.e+"_"+str(num)
				destination_coords_dump  = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/coords_"+str(num)
				destination_lattice_dump = str(U)+"/DOP_"+str(args.dop)+"/"+str(temp)+"/lattice_dump_"+str(num)

				if source_energy_dump == destination_energy_dump or source_coords_dump == destination_coords_dump or source_lattice_dump == destination_lattice_dump:
					continue 
				else:
					try:
				            shutil.copyfile(source_energy_dump,  destination_energy_dump  )
				            shutil.copyfile(source_coords_dump,  destination_coords_dump  )
				            shutil.copyfile(source_lattice_dump, destination_lattice_dump )
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



