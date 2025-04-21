#!/home/satyend/.conda/envs/phase/bin/python

import sys
sys.path.append("/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Analytical_Solvation/exec/new/build")
import numpy as np
import partition_module as pm
import mpmath
import time
import pickle
np.set_printoptions(threshold=10)

if __name__=="__main__":

	print(f"Set up the object.", flush=True)
	partition = pm.Partition(
		Nm=40, 
		Nmm_max=273,
		T=275,
		pw=0.1846742737791685,
		EMM_A=1.3725828290362372,
		EMM_N=1.3725828290362372,
		EMS_A=-9.439946382902741,
		EMS_N=-8.89602248425938,
		alpha=7.333501277066573,
		beta=34.26095605400542,
		dfile="thermo.dump"
	)
	opening_tiles = f'''
	Nm      = {partition.Nm}
	Nmm_max = {partition.Nmm_max}
	T       = {partition.T}
	pw      = {partition.pw}
	EMM_A   = {partition.EMM_A}
	EMM_N   = {partition.EMM_N}
	EMS_A   = {partition.EMS_A}
	EMS_N   = {partition.EMS_N}
	z       = {partition.z}
	ccc     = {partition.contact_constraint_constant}
	alpha   = {partition.alpha}
	beta    = {partition.beta}
	loc     = {partition.loc}
	scale   = {partition.scale}
	thermo  = {partition.dumpfile}
	'''
	f = open("test.dump", 'w')
	print(opening_tiles)
	T_arr = np.linspace(275, 360, 40)
	for T in T_arr:
		partition.T = T
		partition.get_partition_weights_quick()
		partition.get_partition()
		partition.get_free_energy()
		partition.get_average_energy()
		partition.get_average_Nmm()
		partition.get_average_Nmm_a()
		partition.get_average_Nms_a()
		partition.get_Cv()
		f.write(f"{T} | {partition.ave_Nmm} | {partition.ave_Nms_a}\n")
	f.close()

