#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Analyze polymer for rouse modes.")
parser.add_argument ("-N", dest='N', type=int, action='store', help="Length of polymer.")
parser.add_argument ("-p", dest='p', type=int, action='store', help="Number of rouse mode.")
args = parser.parse_args ()

def f_exact (n_, t_, k_, b_, p_, x_):

	t1 = 2 * np.log (x_) 
	t2 = - 6 * (np.sin(p_*np.pi/(2*n_))**2) / (b_ ** 2) * (x_ ** 2)
	t3 = np.log(24 * np.sqrt(6/np.pi) * (np.sin (np.pi*p_/(2*n_))/b_)  ** 3)

	return -k_*t_*(t1+t2+t3)


def Xp (traj, p, N):
	xp = []
	i_list = np.arange (0, N)
	print (i_list + 1/2)
	cos_postfix = np.cos (p*np.pi/N * (i_list+1/2))
	print (f"N = {N}")
	prefactor = np.sqrt(2/N)

	keys_to_probe = list (traj.keys())[1000:]

	for key in keys_to_probe:
		if len(traj[key])!=21:
			print (f"traj[key].shape = {traj[key].shape}")
			print ("Something's fucked.")
			exit ()
		# x1_loop = 0
		# for j in range(N):
		# 	x1_loop += prefactor * traj[key][j,:] * np.cos (p*np.pi/N * (j+1/2))
		# x1_l_mag = np.linalg.norm (x1_loop)

		x1_vec = np.sum (prefactor * traj[key] * cos_postfix[:, np.newaxis], axis=0)
		
		x1_v_mag = np.linalg.norm (x1_vec)

		# print (f"x1_vec = {x1_vec}, x1_v_mag = {x1_v_mag}")
		# if abs(x1_l_mag - x1_v_mag) > 1e-6:
		# 	print ("Somethings off.")
		# 	print (f"Error = {abs(x1_loop-x1_vec)}.")
		# 	exit  ()

		xp.append ( np.linalg.norm (x1_v_mag) )

	return xp



if __name__=="__main__":

	f = open ("coords.nvt.lammpstrj", 'r')
	coords_flag    = False
	timestep_flag  = False
	trajectory = {}

	DOP = args.N
	p   = args.p

	for line in f:
		
		if re.findall ("^ITEM: TIMESTEP$", line):
			polymer = np.zeros ((DOP,3))
			timestep_flag     = True
			coords_flag       = False
			continue
		
		elif re.findall ("^[0-9]+$", line) and timestep_flag:
			timestep = int (line)
			continue
		
		elif re.findall ("^ITEM: NUMBER OF ATOMS$", line):
			timestep_flag = False
			continue

		elif re.findall ("^ITEM: ATOMS id type x y z$", line):
			coords_flag   = True
			continue
		
		elif re.findall ("^.*\\s.*\\s.*\\s.*\\s.*$", line) and coords_flag:
			info = line.strip().split()
			polymer[int(info[0])-1, :] = np.array ([float(info[2]), float(info[3]), float(info[4])])
			if int(info[0]) == DOP:
				trajectory [timestep] = polymer
			continue
		
		else:
			continue


	# get rouse modes 
	# get a list of R coordinates i.e. polymer 
	# given a p value, get the cosine values 
	x1 = Xp (trajectory, p, DOP)
	fig1 = plt.figure()
	ax1  = plt.axes  ()
	(n_h, bins_h, patches) = ax1.hist (x1, density=True, bins=1000)
	fig1.savefig ("x1_hist", dpi=1200, bbox_inches="tight")

	print(f"len(n_h) = {len(n_h)}")
	print(f"len(bins_h) = {len(bins_h)}")

	fig2 = plt.figure(figsize=(5,5))
	ax2  = plt.axes  ()
	F_x1 = -2/3*np.log(n_h)

	b = 1
	T = 2/3
	k = 1
	x_id   = np.linspace (0.001, 20, 1000)
	F_id   = f_exact (DOP, T, k, b, p, x_id)
	F_0  = np.min (F_id)
	F_id   = F_id - F_0

	ax2.plot (x_id, F_id, markersize=0, lw=2, label=f"exact")

	F_x1 = F_x1 - np.min (F_x1)

	ax2.plot ((bins_h[0:-1]+bins_h[1:])/2, F_x1, label=f"unbiased")
	ax2.legend(loc="upper right")
	ax2.set_xlim (0, 20)
	ax2.set_ylim (0, 6)
	ax2.set_xticks (np.linspace (0, 20, 5))
	ax2.set_yticks (np.linspace (0, 6, 4) )
	fig2.savefig ("f_x1", dpi=1200, bbox_inches="tight")



