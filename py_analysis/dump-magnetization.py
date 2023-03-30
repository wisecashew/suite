#!/home/satyend/.conda/envs/data_analysis/bin/python

import pandas as pd 
import numpy as np 
import re
import aux
import argparse 
import os 

parser = argparse.ArgumentParser (description="Dump magnetization.")
parser.add_argument ("--file-path", dest='fp', action='store', type=str, help="Enter a file path to orientation file.")

args = parser.parse_args()



def dump_magnetization ( ortn_file_path, starting_index ):
	Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
            3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
            6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
            9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
            12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
            15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
            18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
            24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }

	start_str = "Dumping coordinates at step"
	end_str   = "END"
	magnetization_equi  = []
	magnetization = []
	num_list = aux.dir2nsim (os.listdir(ortn_file_path))
	for idx in num_list:
		f = open(ortn_file_path+"/coords_"+str(idx)+".mc", 'r')

		extract_orr = False
		start_bool  = False
		monomer_order  = np.array([0,0,0])
		magnetization.clear()

		for line in f:
			if re.match ( start_str, line ):
				a = re.search ("\d+", line)
				# extract_orr = True
				monomer_order = np.array([0,0,0], dtype=np.float64)
				if int ( a.group(0) ) == 0:
					start_bool = True
				count          = 0

			elif re.match ("START", line):
				extract_orr = True

			elif re.match ( end_str, line ) and start_bool:
				magnetization.append (np.linalg.norm(monomer_order))
				monomer_order = np.array ([0,0,0])
				extract_orr = False

			elif extract_orr and start_bool:
				# print (line, )
				# print (extract_loc_from_string(line))
				monomer_or = aux.extract_loc_from_string ( line ) [3]
				count = 0
				monomer_order += Or2Dir[monomer_or]
		magnetization_equi.extend (magnetization[-starting_index:])
		f.close()
	d = {"step": [1], "magnetization": [np.mean(magnetization_equi)] }
	df = pd.DataFrame.from_dict(d)
	df.to_csv (ortn_file_path+"/magnetization_dump.mc", index=False, sep='|')

	return



if __name__=="__main__":

	dump_magnetization (args.fp, 2000)

