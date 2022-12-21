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
parser.add_argument('-pdb', dest='pdb', metavar='a.pdb', action='store', type=str, help='Enter a name for the pdb file.')

args = parser.parse_args() 

def location (lattice_index, x, y, z):
    zc = lattice_index // (z*z)
    yc = (lattice_index % (z*z)) // y
    xc = ( ( lattice_index % (z*z) ) % y ) % x
    return np.array([xc, yc, zc])

if __name__ == "__main__":

	start = time.time()
	step  = args.s
	start_str = "FINAL STEP: " + str(args.s)
	print (start_str)
	end_str   = "END."
	x = args.x
	y = args.y
	z = args.z
	frac = aux.get_frac(args.ff)
	nsol2 = int(np.floor(((x**3)-args.dop)*frac))
	nsol1 = int(x**3 - args.dop -nsol2)
	print ("nsol1 = ",nsol1, ", nsol2 =", nsol2)
	step_bool = False

	f = open (args.lat, 'r')
	for line in f:
		if re.findall (start_str, line):
			r    = re.findall ("\d+", line) 
			step = r[0]
			step_bool = True
			polymer = np.zeros((args.dop, 3))
			solvent1 = np.zeros((nsol1,   3))
			solvent2 = np.zeros((nsol2,   3))
			m_num  = 0
			s1_num = 0
			s2_num = 0
		elif re.findall (end_str, line) and step_bool:
			break
		elif step_bool:
			info = line.strip().split() 
			if info[1] == "m1,":
				# polymer[m_num] = location   (int(info[2]), x, y, z)
				# m_num  += 1
				continue
			elif info[1] == "s1,":
				solvent1[s1_num] = location (int(info[2]), x, y, z)
				s1_num += 1
			elif info[1] == "s2,":
				solvent2[s2_num] = location (int(info[2]), x, y, z)
				s2_num += 1
	if not step_bool:
		print ("step not found. exiting...", flush=True)
		print (step_bool)
		exit()
	
	f.close()
	######################################
	pdict = aux.get_pdict (args.c, int(step), 32, x, y, z)
	polymer = pdict[int(step)][0]
	polymer = aux.unfuck_polymer(polymer, x, y, z)


	#######################################
	# polymer = aux.unfuck_polymer (polymer, x, y, z)
	# com = np.mean(polymer, axis=0)
	# print ("com = ", com)
	# polymer = polymer-com
	# solvent1 = solvent1-com
	# solvent2 = solvent2-com
	f = open (args.pdb, 'w')

	f.write("COMPND    MY_POLYMER\n")
	f.write("AUTHOR    SAT\n")
	count = 1 
	for row in polymer:
		f.write ("ATOM  {:>5} {:>2}{:<2} POL P{:>4} {:>8}{:>8}{:>8}{:>6}{:>6}      {:>4}{:<1}\n".format( count, "M", " ", 1, row[0]*10, row[1]*10, row[2]*10, 1.00, 1.00, "M", "M" ) )
		count+=1

	j = 2
	scount = copy.copy(count)
	for sol in solvent1:
		f.write ("ATOM  {:>5} {:>2}{:<2} SOL  {:>4} {:>8}{:>8}{:>8}{:>6}{:>6}      {:>4}{:<1}\n".format( scount, "S", "1", j, sol[0]*10, sol[1]*10, sol[2]*10, 1.00, 1.00, "S", "1" ) )
		scount += 1
	j += 1
	
	for sol in solvent2:
		f.write ("ATOM  {:>5} {:>2}{:<2} SOL  {:>4} {:>8}{:>8}{:>8}{:>6}{:>6}      {:>4}{:<1}\n".format( scount, "S", "2", j, sol[0]*10, sol[1]*10, sol[2]*10, 1.00, 1.00, "S", "2" ) )
		scount += 1
	j += 1

	for i in range(1, count):
		if i == 1:
			f.write ("CONECT{:>5}{:>5}\n".format(i, i+1))
		elif i == count-1:
			f.write ("CONECT{:>5}{:>5}\n".format(i, i-1))
		else:
			f.write("CONECT{:>5}{:>5}{:>5}\n".format(i, i-1, i+1 ))

	f.close()
	stop = time.time()
	print ("Run time for N = " + str(args.dop) + " is {:.2f} seconds.".format(stop-start), flush=True)

