#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import nglview as nv
import argparse
import pickle
import MDAnalysis.transformations as xforms
from mpl_toolkits.mplot3d import Axes3D
import latt_lmp as ll

parser = argparse.ArgumentParser(description="Compute hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--hpkl",         dest="hpkl", type=str, action="store", help="enter address of hydrogen bonding pickle dump.", required=True)
parser.add_argument("--Lpkl",         dest="lpkl", type=str, action="store", help="enter address of final processed lattice.",      required=True)
args = parser.parse_args()

if __name__=="__main__":

	# start the clock
	start = time.time()

	# get the lattice
	my_lattice = open(args.lpkl, "rb")
	lattice    = pickle.load(my_lattice)
	my_lattice.close()

	# get the number of aligned water solvents
	for idx in lattice.lattice:
		# indices of the oxygen and nitrogen engaging in hbonding
		rel_idx_O = lattice.hba['Op_with_water'][0][:,0] == idx
		rel_idx_N = lattice.hba['Np_with_water'][0][:,0] == idx

		# get the water indices
		waterO_wOp = lattice.hba["Op_with_water"][0][rel_idx_O][:,1]
		waterO_wNp = lattice.hba["Np_with_water"][0][rel_idx_N][:,1]

		common_elements = np.intersect1d(waterO_wOp, waterO_wNp)
		# print(f"{idx}. common elements = {common_elements}")

		for ce in common_elements:
			Hidx_O = lattice.hba["Op_with_water"][0][rel_idx_O][lattice.hba["Op_with_water"][0][rel_idx_O][:,1] == ce][:,2]
			Hidx_N = lattice.hba["Np_with_water"][0][rel_idx_N][lattice.hba["Np_with_water"][0][rel_idx_N][:,1] == ce][:,2]
			if Hidx_O == Hidx_N:
				print(f"Same H bonding with polymer O and N!")

		Nmsa = len(lattice.lattice[idx]['solvent'])
		Ns   = np.sum(lattice.hba['Op_with_water'][0][:,0] == idx)+np.sum(lattice.hba['Np_with_water'][0][:,0] == idx)
		if Nmsa != Ns:
			print(f"Mismatch on idx {idx}!")
		# print(f"Number of aligned solvent particles is {Nmsa}", end=' | ')
		# print(f"Counting participating water molecules = {Ns}")


	# plot the first timestep on the lattice
	fig = plt.figure()
	ax  = fig.add_subplot(111, projection="3d")

	for solvent in lattice.lattice[0]["solvent"]:
		ll.plot_cube(ax, solvent, 1, "steelblue")

	for monomer in lattice.lattice[0]["polymer"]:
		ll.plot_cube(ax, monomer, 1, "coral")

	ax.set_xlim(-20, 20)
	ax.set_ylim(-20, 20)
	ax.set_zlim(-20, 20)

	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False

	ax.grid(False)

	fig.savefig("ts0.png", dpi=1200, bbox_inches="tight")
	"""