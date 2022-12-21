#!/usr/licensed/anaconda3/2020.7/bin/python

import re 
import numpy as np 
import argparse 
import aux
import os

parser = argparse.ArgumentParser(description="Read a trajectory file and obtain a radius of gyration plot given a degree of polymerization over a range of temperatures and potential energy surfaces.")
parser.add_argument('-s', dest='s', metavar='filename', action='store', type=str, help='Name of single run.')
parser.add_argument('-d', dest='d', metavar='filename', action='store', type=str, help='Name of direct run.')
parser.add_argument('-m', dest='m', metavar='filename', action='store', type=str, help='Name of merged run.')
args = parser.parse_args()



SMALL_RUN_FILE = open (args.s, 'r') 
SMALL_RUN_DICT = {} 
U_re     = "U ="
flory_re = "flory exp:"
T_re     = "T:"
r2_re    = "r2:" 

U_master_list = aux.dir2U ( os.listdir (".") )

U_list = [] 

for line in SMALL_RUN_FILE:
	if re.findall (U_re, line):
		U_id = re.findall("U\d+", line) 
		U_list.append ( U_id[0] )
		SMALL_RUN_DICT[U_id[0]] = {} 
		continue
	elif re.findall (flory_re, line):
		flory_id = line.strip().split()[2:] # re.findall ("\d+\.\d+", line)
		continue
	elif re.findall (r2_re, line):
		r2_id    = line.strip().split()[1:] # re.findall ("\d+\.\d+", line)
		continue
	elif re.findall (T_re, line):
		T_id    = line.strip().split()[1:] # re.findall ("\d+\.\d+", line)
		T_list  = T_id 
		for i in range (len(T_list)):
			SMALL_RUN_DICT[U_id[0]][T_list[i]] = (float(flory_id[i]), float(r2_id[i]))

SMALL_RUN_FILE.close() 

LARGE_RUN_FILE = open (args.d, 'r')
LARGE_RUN_DICT = {} 

U_list = [] 
T_list = [] 

for line in LARGE_RUN_FILE:
	if re.findall (U_re, line):
		U_id = re.findall ("U\d+", line) 
		U_list.append (U_id[0])
		continue
	elif re.findall (flory_re, line):
		flory_id = line.strip().split()[2:]
		continue
	elif re.findall (r2_re, line):
		r2_id   = line.strip().split()[1:] 
		continue 
	elif re.findall (T_re, line):
		T_id   = line.strip().split()[1:] 
		T_list = T_id 
		LARGE_RUN_DICT[U_id[0]] = {} 
		for i in range (len(T_list)):
			LARGE_RUN_DICT[U_id[0]][T_list[i]] = (float(flory_id[i]), float(r2_id[i])) 

LARGE_RUN_FILE.close()

# start rewriting file 

if ( U_list != U_master_list ):
	print ("Something is not right with the set-up.")
	exit ()

new_file = open (args.m, 'w') 

for U in U_master_list:
	new_file.write ("U = " + U + ":\n")
	new_file.write ("flory exp:  ")
	for T in T_list:

		if (not U in SMALL_RUN_DICT):
			new_file.write (str(LARGE_RUN_DICT[U][T][0])+"  ")
		elif (not T in SMALL_RUN_DICT[U] ):
			new_file.write (str(LARGE_RUN_DICT[U][T][0])+"  ")
		elif ( float(SMALL_RUN_DICT[U][T][1]) >= float(LARGE_RUN_DICT[U][T][1]) ):
			new_file.write (str(SMALL_RUN_DICT[U][T][0])+"  ")
		elif ( float(SMALL_RUN_DICT[U][T][1]) < float(LARGE_RUN_DICT[U][T][1]) ):
			new_file.write (str(LARGE_RUN_DICT[U][T][0])+"  ")

	new_file.write ("\n")
	new_file.write ("r2:         ")
	for T in T_list:
		if ( not U in SMALL_RUN_DICT ):
			new_file.write (str(LARGE_RUN_DICT[U][T][1])+"  ")
		elif ( not T in SMALL_RUN_DICT[U] ):
			new_file.write (str(LARGE_RUN_DICT[U][T][1])+"  ")
		elif ( float(SMALL_RUN_DICT[U][T][1]) >= float(LARGE_RUN_DICT[U][T][1]) ):
			new_file.write (str(SMALL_RUN_DICT[U][T][1])+"  ")
		elif ( float(SMALL_RUN_DICT[U][T][1]) < float(LARGE_RUN_DICT[U][T][1]) ):
			new_file.write (str(LARGE_RUN_DICT[U][T][1])+"  ")
	new_file.write ("\n")
	new_file.write ("T:          ")
	for T in T_list:
		new_file.write ((T + "  "))
	new_file.write ("\n")

new_file.close()

