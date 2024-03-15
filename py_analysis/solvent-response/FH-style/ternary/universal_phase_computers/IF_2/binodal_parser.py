import numpy as np
import argparse
import time
import warnings
import linecache
import ternary
import pickle
import os

import argparse
parser = argparse.ArgumentParser(description='Create a skeleton solution for the binodal. This is a memory-intensive computation.'         )
parser.add_argument('--addr-get-from', metavar='bpkl',      dest='aget',     type=str, action='store', help="Name of binodal pkl directory.", default=None )
parser.add_argument('--no-rtw',         dest='nrtw',       action='store_true',  help="Dont print out the runtime warning.", default=False)
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

#########################################

##########################################

if __name__=="__main__":

	start = time.time()

	#=================================
	l_binodals = os.listdir(args.aget)
	print(f"Number of binodals is {len(l_binodals)}.", flush=True)

	for bfyle in l_binodals:
		print(f"Probing {args.aget}/{bfyle}...")
		f = open(args.aget+"/"+bfyle, 'rb')
		BINODALS = pickle.load(f)
		f.close()

		for idx, key in enumerate(BINODALS["hull_info"]["binodal"]):
			B1 = key[0][0]
			B2 = key[0][1]
			B1, keep = ternary.remove_close_rows(B1, 1e-6)
			B2       = B2[keep]

			print(f"len(B1) = {len(B1)}, len(B2) = {len(B2)}", flush=True)

			c1  = BINODALS["groupings"][(0,0)]["raw_crits"][0]
			c2  = BINODALS["groupings"][(0,0)]["raw_crits"][1]
			mid = (c1+c2)/2
			mid_axis    = (mid - c1)[0:2]/np.linalg.norm((mid-c1)[0:2])
			adj_neg_arm = (B1[:,0:2]-mid[0:2])/np.linalg.norm(B1[:,0:2]-mid[0:2], axis=1)[:,np.newaxis]
			angles      = np.arccos(np.sum(adj_neg_arm * mid_axis, axis=1))
			B1  = B1[np.argsort(angles)]
			B2  = B2[np.argsort(angles)]

			BINODALS["hull_info"]["binodal"][idx][0][0] = B1
			BINODALS["hull_info"]["binodal"][idx][0][1] = B2

			print(f"BINODALS[\"hull_info\"][\"binodal\"][{idx}][0][0] = {len(BINODALS['hull_info']['binodal'][idx][0][0])}")
			print(f"BINODALS[\"hull_info\"][\"binodal\"][{idx}][0][1] = {len(BINODALS['hull_info']['binodal'][idx][0][1])}")
		
		print(f"Probed! Dumping now.", flush=True)
		f = open(args.aget+"/"+bfyle[:-4]+".mod.pkl", 'wb')
		pickle.dump(BINODALS, f)
		f.close()

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

