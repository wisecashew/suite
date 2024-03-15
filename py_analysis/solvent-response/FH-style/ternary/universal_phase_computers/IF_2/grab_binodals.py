import time
import argparse
import warnings
import linecache
import random
import os
import re
import shutil

import argparse
parser = argparse.ArgumentParser(description='This will make images of binodals.')
parser.add_argument('--addr-get-from', metavar='bpkl',      dest='aget',     type=str, action='store', help="Name of binodal pkl directory.", default=None )
parser.add_argument('--addr-send-to-chips',  metavar='bpkl',      dest='asendps',    type=str, action='store', help="Name of binodal pkl directory to dump to.", default=None )
parser.add_argument('--addr-send-to-chipc',  metavar='bpkl',      dest='asendpc',    type=str, action='store', help="Name of binodal pkl directory to dump to.", default=None )
parser.add_argument('--addr-send-to-chisc',  metavar='bpkl',      dest='asendsc',    type=str, action='store', help="Name of binodal pkl directory to dump to.", default=None )
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

if __name__=="__main__":

	start = time.time()

	#####################################################
	l_binodals   = os.listdir(args.aget)
	print(f"Number of binodals is {len(l_binodals)}.", flush=True)

	chips_emph = []
	chipc_emph = []
	chisc_emph = []
	for num, file in enumerate(l_binodals):

		if re.findall("mod", file):
			# get the parameters
			chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1))
			chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1))
			chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1))
			vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1))
			vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1))
			vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1))
			if chips < -1:
				chips_emph.append(file)
			elif chipc < -1:
				chipc_emph.append(file)
			elif chisc < -1:
				chisc_emph.append(file)
		else:
			continue

	# choose 100 binodals form each section
	n = 100
	chips_to_grad = [chips_emph[i] for i in random.sample(range(len(chips_emph)), n)]
	chipc_to_grad = [chipc_emph[i] for i in random.sample(range(len(chipc_emph)), n)]
	chisc_to_grad = [chisc_emph[i] for i in random.sample(range(len(chisc_emph)), n)]

	# copy files from a binodals 
	for cps in chips_to_grad:
		shutil.copyfile(f"{args.aget[:-1]}{cps}", f"{args.asendps[:-1]}{cps}")
	
	for cpc in chipc_to_grad:
		shutil.copyfile(f"{args.aget[:-1]}{cpc}", f"{args.asendpc[:-1]}{cpc}")

	for csc in chisc_to_grad:
		shutil.copyfile(f"{args.aget[:-1]}{csc}", f"{args.asendsc[:-1]}{csc}")

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

