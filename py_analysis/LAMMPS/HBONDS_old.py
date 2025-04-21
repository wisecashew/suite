#!/home/satyend/.conda/envs/phase/bin/python

import lmp_settings
import lmp_data
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import argparse

parser = argparse.ArgumentParser(description="Compute hydrogen bonding information specifically for REMD SPCE+NIPAM FILES!")
parser.add_argument("--datafile",     dest="df",   type=str, action="store", help="enter address of datafile.",        required=True)
parser.add_argument("--settingsfile", dest="sf",   type=str, action="store", help="enter address of settingsfile.",    required=True)
parser.add_argument("--coords",       dest='c',    type=str, action="store", help="enter address of coordinate file.", required=True)
parser.add_argument("--pkl",          dest="pkl",  type=str, action="store", help="enter address of pickle dump.",     required=True)
args = parser.parse_args()

if __name__=="__main__":

	start = time.time()

	# get the names of the files
	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	# generate the universe object
	print(f"Making the universe...", end=' ', flush=True)
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"made!", flush=True)

	H_polymer = "type 7"
	H_water   = "type 11"
	O_polymer = "type 5"
	N_polymer = "type 6"
	O_water   = "type 10"

	# access atom types for hbonding
	print(f"Run hbond analysis...", flush=True)
	hbonds_with_Opoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=O_polymer, donors_sel=O_water)
	hbonds_with_Npoly  = HBA(universe=u, hydrogens_sel=H_water,   acceptors_sel=N_polymer, donors_sel=O_water)
	hbonds_with_Owater = HBA(universe=u, hydrogens_sel=H_polymer, acceptors_sel=O_water,   donors_sel=N_polymer)
	hbonds_with_Op_Np  = HBA(universe=u, hydrogens_sel=H_polymer, acceptors_sel=O_polymer, donors_sel=N_polymer)

	# run hbond analysis with MDA
	print(f"Running hbond analysis with water H being donated by water O and accepted by polymer O.", end=' ', flush=True)
	hbonds_with_Opoly.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with water H being donated by water O and accepted by polymer N.", end=' ', flush=True)
	hbonds_with_Npoly.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with amidic H being donated by amidic N and accepted by water O.", end=' ', flush=True)
	hbonds_with_Owater.run()
	print(f"done!", flush=True)
	print(f"Running hbond analysis with amidic H being donated by amidic N and accepted by polymer O.", end=' ', flush=True)
	hbonds_with_Op_Np.run()
	print(f"done!", flush=True)

	# set up the hbonds
	hbonds = [hbonds_with_Opoly, hbonds_with_Npoly, hbonds_with_Owater, hbonds_with_Op_Np]
	with open(args.pkl, "wb") as f:
		pickle.dump(hbonds, f)

	# print out stuff
	print(f"Done running.", flush=True)
	print(f"Amidic N with Water H: \n{hbonds_with_Npoly.results.hbonds.shape}",    flush=True)
	print(f"Carbonyl O with Water H: \n{hbonds_with_Opoly.results.hbonds.shape}",  flush=True)
	print(f"Water O with Amide H: \n{hbonds_with_Owater.results.hbonds.shape}",    flush=True)
	print(f"Polymer H with Polymer O: \n{hbonds_with_Op_Np.results.hbonds.shape}", flush=True)
	print(f"Completed computation.", flush=True)

	# get the times in
	stop = time.time()
	print(f"Computation took {stop-start} seconds.", flush=True)
