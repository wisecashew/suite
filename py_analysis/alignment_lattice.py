#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np 
import re 
import pandas as pd
import aux 
import time 
import sys 
import itertools
import copy

sys.stdout.flush() 

import argparse 
parser = argparse.ArgumentParser(description="Read a trajectory file and obtain the pdb file for the polymer.")
parser.add_argument('-dop', metavar='DOP', dest='dop', type=int, action='store', help='enter a degree of polymerization.')
parser.add_argument('-s', metavar='S', type=int, dest='s', action='store', help='start parsing after this index.', default=100)
parser.add_argument('--lat', dest='lat', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('-c', dest='c', metavar='coords.txt', action='store', type=str, help='Name of energy dump file to parse information.', default='coords.txt')
parser.add_argument('-ff', dest='ff', action='store', type=str, help='Address to geom_and_esurf.txt')
parser.add_argument('-x', dest='x', metavar='X', action='store', type=float, help='Enter x-dimension of box.')
parser.add_argument('-y', dest='y', metavar='Y', action='store', type=float, help='Enter y-dimension of box.')
parser.add_argument('-z', dest='z', metavar='Z', action='store', type=float, help='Enter z-dimension of box.')

args = parser.parse_args() 

Or2Dir = { 0: np.asarray([1,0,0]), 1: np.asarray ([0,1,0]), 2: np.asarray([0,0,1]), \
			3: np.asarray([-1,0,0]), 4: np.asarray([0,-1,0]), 5: np.asarray([0,0,-1]), \
			6: np.asarray([1/np.sqrt(2),1/np.sqrt(2),0]), 7: np.asarray([1/np.sqrt(2), 0, 1/np.sqrt(2)]), 8: np.asarray ([1/np.sqrt(2),-1/np.sqrt(2),0]), \
			9: np.asarray([1/np.sqrt(2),0,-1/np.sqrt(2)]), 10: np.asarray([-1/np.sqrt(2),1/np.sqrt(2),0]), 11: np.asarray([-1/np.sqrt(2),0,1/np.sqrt(2)]), \
			12: np.asarray([-1/np.sqrt(2),-1/np.sqrt(2),0]), 13: np.asarray([-1/np.sqrt(2),0,-1/np.sqrt(2)]), 14: np.asarray([0,1/np.sqrt(2),1/np.sqrt(2)]), \
			15: np.asarray([0,1/np.sqrt(2),-1/np.sqrt(2)]), 16: np.asarray([0,-1/np.sqrt(2),1/np.sqrt(2)]), 17: np.asarray([0,-1/np.sqrt(2),-1/np.sqrt(2)]), \
			18: np.asarray([1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 19: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 20: np.asarray([1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
			21: np.asarray([1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]), 22: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]), 23: np.asarray([-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)]), \
			24: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)]), 25: np.asarray([-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)]) }

v = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], \
	[0, -1, 0], [0, 0, 1], [0, 0, -1], \
	[1, 1, 0], [-1, 1, 0], [1, -1, 0], \
	[-1, -1, 0], [0, 1, 1], [0, -1, 1], \
	[0, 1, -1], [0, -1, -1], \
	[1, 0, 1], [-1, 0, 1], [1, 0, -1], \
	[-1, 0, -1], [1, 1, 1], \
	[-1, 1, 1], [1, -1, 1], \
	[1, 1, -1], [-1, -1, 1], \
	[-1, 1, -1], [1, -1, -1], [-1, -1, -1]]
v = [np.asarray(u) for u in v] 
v = np.asarray (v) 

def alignment_condition(o1, o2):
	dot_prod = np.dot(Or2Dir[o1], Or2Dir[o2])
	if dot_prod > 0.54:
		return True
	else:
		return False

def lat2location (lattice_index, x, y, z):
    zc = lattice_index // (z*z)
    yc = (lattice_index % (z*z)) // y
    xc = ( ( lattice_index % (z*z) ) % y ) % x
    return np.array([xc, yc, zc])

def location2latidx(location, x, y, z):
	return int(location[2]*z*z + location[1]*y + location[0])

def impose_pbc(loc, x, y, z):
	loc[0] = (( loc[0]%x ) +x ) %x
	loc[1] = (( loc[1]%y ) +y ) %y
	loc[2] = (( loc[2]%z ) +z ) %z
	return loc 

def obtain_ne_list(location, x, y, z):
	ne_list = []
	for i in range(26):
		ne = location + v[i]
		impose_pbc(ne, x, y, z)
		ne = location2latidx(ne, x, y, z)
		ne_list.append(ne)
	return ne_list

if __name__ == "__main__":

	start = time.time()
	step  = args.s
	start_str = "FINAL STEP: " + str(args.s)
	end_str   = "END."
	x = args.x
	y = args.y
	z = args.z
	frac = float(aux.get_frac(args.ff))
	print(type(x), type(args.dop), type(frac))
	nsol2 = int(np.floor(((x**3)-args.dop)*frac))
	nsol1 = int(x**3 - args.dop - nsol2)
	print ("nsol1 = ",nsol1, ", nsol2 =", nsol2)
	step_bool = False

	# create containers
	LATTICE   = dict()
	polymer   = []
	solvent   = []
	cosolvent = []

	f = open (args.lat, 'r')
	for line in f:
		if re.findall (start_str, line):
			r    = re.findall ("\d+", line) 
			step = r[0]
			step_bool = True
			m_num  = 0
			s1_num = 0
			s2_num = 0
			print("line 1: ", line)
		elif re.findall (end_str, line) and step_bool:
			break
		elif step_bool:
			info = line.strip().split(",") 
			# print(info)
			LATTICE[int(info[2].strip())] = (info[1].strip(), int(info[0].strip()))
			# print(info[1])
			if info[1].strip() == "m1":
				# print(f"m1: ", info)
				polymer.append(int(info[2]))
			elif info[1].strip() == "s1":
				# print(f"s1: ", info)
				solvent.append(int(info[2]))
			elif info[1].strip() == "s2":
				# print(f"s2: ", info)
				cosolvent.append(int(info[2]))
	if not step_bool:
		print ("step not found. exiting...", flush=True)
		print (step_bool)
		exit()
	
	f.close()

	# print("polymer: ", polymer)
	# print("solvent: ", solvent)
	# print("cosolvent: ", cosolvent)

	# go to each monomer location and get the orientation
	# then, obtain all the "solvent" (S1) neighbors of the polymer (M1) and identify all the S1 aligned with M1
	# then, go to all the S1 neighbors of the aligned S1 and identify all the S2 neighbors that are aligned with S1
	# return all instances of aligned M1-S1-S2 triplicates
	aligned_triples = []
	for monomer in polymer:
		monomer_loc = lat2location(monomer, x, y, z)
		monomer_or  = LATTICE[monomer][1]
		monomer_ne_list     = obtain_ne_list(monomer_loc, x, y, z)
		for neighbor in monomer_ne_list:
			if neighbor in solvent:
				neighbor_loc = lat2location(neighbor, x, y, z)
				neighbor_or  = LATTICE[neighbor][1]
				if alignment_condition(monomer_or, neighbor_or):
					# get all the neighbors of the solvent
					solvent_ne_list = obtain_ne_list(neighbor_loc, x, y, z)
					for solvent_neighbor in solvent_ne_list:
						if solvent_neighbor in cosolvent:
							solvent_neighbor_loc = lat2location(solvent_neighbor, x, y, z)
							solvent_neighbor_or  = LATTICE[solvent_neighbor][1]
							if alignment_condition(neighbor_or, solvent_neighbor_or):
								aligned_triples.append([monomer, neighbor, solvent_neighbor])
	
	print (aligned_triples)
	#######################################
	# pdict = aux.get_pdict (args.c, int(step), 32, x, y, z)
	# polymer = pdict[int(step)][0]
	# polymer = aux.unfuck_polymer(polymer, x, y, z)


