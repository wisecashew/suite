#!/usr/licensed/anaconda3/2020.7/bin/python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import StrMethodFormatter
import argparse
import copy
import sys
sys.path.insert (0, '/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/py_analysis')
import aux
import os
import re

parser = argparse.ArgumentParser(description="Supply a CG FF.")
parser.add_argument ("--form", dest='form', type=int, help="Select CG FF form.")
args = parser.parse_args()

if __name__=="__main__":
	x = y = z = 34
	T_list = aux.dir2float (os.listdir("."))
	T_list = [0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 1.0]
	print (T_list)
	E_ms = 0
	if args.form == 2:
		for T in T_list:
			kT = float (T)
			try:
				f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
				if kT >= 20.0:
					E_mm_a = -0.5
					E_mm_n = -0.25
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_a = {E_mm_a}\n")
					f.write (f"Emm_n = {E_mm_n}\n")
					f.write (f"Ems   = 0\n")
					f.write ("END OF FILE")
				elif kT >= 0.5:
					E_mm_a = -0.5
					E_mm_n = -0.25
					f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_a = {E_mm_a}\n")
					f.write (f"Emm_n = {E_mm_n}\n")
					f.write (f"Ems   = 0\n")
					f.write ("END OF FILE")
				else:
					E_mm_a = 0.25 * kT
					E_mm_n = 0.5 * kT
					f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_a = {E_mm_a}\n")
					f.write (f"Emm_n = {E_mm_n}\n")
					f.write (f"Ems   = 0\n")
					f.write ("END OF FILE")
				f.close ()
			except FileNotFoundError:
				continue
	elif args.form == 3:
		for T in T_list:
			kT = float (T)
			try:
				f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
				if kT >= 20.0:
					E_mm_1 = -0.5
					E_mm_2 = -0
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_1 = {E_mm_1}\n")
					f.write (f"Emm_2 = {E_mm_2}\n")
					f.write ("END OF FILE")
				elif kT >= 0.5:
					E_mm_1 = -0.5
					E_mm_2 = 0
					f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_1 = {E_mm_1}\n")
					f.write (f"Emm_2 = {E_mm_2}\n")
					f.write ("END OF FILE")
				else:
					E_mm_1 = 0.5 * kT
					E_mm_2 = 0.25 * kT
					f = open (str(T) + "/FORM" + str(args.form) + "/MODEL1/geom_and_esurf.txt", 'w')
					f.write (f"x = {x}\n")
					f.write (f"y = {y}\n")
					f.write (f"z = {z}\n")
					f.write (f"kT = {kT}\n")
					f.write (f"Emm_1 = {E_mm_1}\n")
					f.write (f"Emm_2 = {E_mm_2}\n")
					f.write ("END OF FILE")
				f.close ()
			except FileNotFoundError:
				continue
	else:
		print ("Invalid form provided: "+str(args.form))
