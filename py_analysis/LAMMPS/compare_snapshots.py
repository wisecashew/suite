#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
import numpy as np
import argparse
import time

parser = argparse.ArgumentParser(description="Compare configurations.")
parser.add_argument("--equal",        dest="eq",    type=str, action="store", help="enter equal tuples.")
parser.add_argument("--datafile",     dest="df",    type=str, action='store', help="enter address of datafile.")
parser.add_argument("--settingsfile", dest="sf",    type=str, action='store', help="enter address of settingsfile.")
parser.add_argument("--coords",       dest='c',     type=str, action='store', help="enter address of coordinate file.")
args = parser.parse_args()

def equality_frames(frame1, frame2):

	# compare cell dimensions
	if not np.allclose(frame1[0], frame2[0]):
		# print(f"Dimensions don't match.", end=' ', flush=True)
		return False

	else:
		# compare atomic coordinates
		if not np.allclose(frame1[1], frame2[1]):
			# print(f"Atomic coordinates don't match.", end=' ', flush=True)
			return False
		else:
			return True
	return

if __name__=="__main__":

	# get the start time
	start = time.time()

	# set up the inputs
	settings = args.sf
	topo     = args.df
	traj     = args.c
	dt       = 2

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)

	# compare two frames
	ltraj = u.trajectory.n_frames

	# container of trajectories
	frame_container = []
	equal_frames    = []

	print(f"Making the frame_container...", flush=True)

	# store all the info
	for ts in u.trajectory:
		# print(f"@ ts = {ts}", flush=True)
		frame_container.append((np.copy(ts.dimensions), np.copy(ts.positions)))

	print(f"Made the frame_container!", flush=True)
	print(f"Number of frames is {len(frame_container)}.", flush=True)

	# compare configs/frames
	for i in range(0, ltraj):
		for j in range(i+1, ltraj):
			if equality_frames(frame_container[i], frame_container[j]):
				print(f"Comparing frame {i*10000} and {j*10000} are equal!", flush=True)
				equal_frames.append((i, j))
			else:
				pass
	
	# print out frames
	tuple_dump = open(args.eq, 'w')
	for ef in equal_frames:
		tuple_dump.write(f"{ef[0]*10000} {ef[1]*10000}\n")
	tuple_dump.close()

	# get the end time
	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)
