#!/home/satyend/.conda/envs/ML/bin/python

import numpy as np
import argparse
import time
import warnings
import linecache
import pickle
import os
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pandas as pd
import mpltern
import ternary
import phase
import re
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

parser = argparse.ArgumentParser(description='Make databases from .pkl files for binodals.')
parser.add_argument('--mesh-density',     metavar='mdense',  dest='meshd',     type=int,   action='store', help='enter mesh density.', default=50)
parser.add_argument('--dirs-of-binodals', metavar='dirs',    dest='dirs',      type=str, action='store',   help='address of binodals. Enter address as follows: \"/path/to/directory\"', nargs='+', default=None)
parser.add_argument('--database-single',  metavar='SDB',     dest='db_single', type=str, action='store',   help='name of database for homoegenous points.')
parser.add_argument('--database-double',  metavar='DDB',     dest='db_double', type=str, action='store',   help='name of database for splitting ponints.' )
parser.add_argument('--database-triple',  metavar='TDB',     dest='db_triple', type=str, action='store',   help='name of database for triple points.'     )
parser.add_argument('--database-comb',    metavar='CDB',     dest='db_comb',   type=str, action='store',   help='name of database for all points.'        )
parser.add_argument('--only-labels',  dest='ol',      action='store_true', help='Make a database with only labels.',       default=False)
parser.add_argument('--with-weights', dest='ww',      action='store_true', help='Make a database with splitting weights.', default=False)
parser.add_argument('--combine',      dest='combine', action='store_true', help='Combine databases.',                      default=False)
parser.add_argument('--nrtw',         dest='nrtw',    action='store_true', help='run time warning display.',               default=False)
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

##########################################

def generate_points_inside_triangle(vertices, N):
    u = np.random.rand(N, 1)
    v = np.random.rand(N, 1)

    is_inside_triangle = u + v <= 1
    u_inside = u[is_inside_triangle]
    v_inside = v[is_inside_triangle]

    w_inside = 1 - u_inside - v_inside

    x = (u_inside * vertices[0, 0] + v_inside * vertices[1, 0] + w_inside * vertices[2, 0]).reshape(-1, 1)
    y = (u_inside * vertices[0, 1] + v_inside * vertices[1, 1] + w_inside * vertices[2, 1]).reshape(-1, 1)

    points_inside_triangle = np.hstack((x, y))

    return points_inside_triangle

##########################################

def double_split(C1, C2, point):
	V1 = (C2[1]*point[0] - C2[0]*point[1])/(C2[1]*C1[0] - C2[0]*C1[1])
	V2 = (C1[1]*point[0] - C1[0]*point[1])/(C2[0]*C1[1] - C2[1]*C1[0])
	V  = V1+V2
	return np.array([V1/V, V2/V]).T

##########################################

def triple_split(C1, C2, C3, point):
	D  = np.array([[C1[0], C2[0], C3[0]], [C1[1], C2[1], C3[1]], [C1[2], C2[2], C3[2]]])
	D1 = np.array([[point[0], C2[0], C3[0]], [point[1], C2[1], C3[1]], [point[2], C2[2], C3[2]]])
	D2 = np.array([[C1[0], point[0], C3[0]], [C1[1], point[1], C3[1]], [C1[2], point[2], C3[2]]])
	D3 = np.array([[C1[0], C2[0], point[0]], [C1[1], C2[1], point[1]], [C1[2], C2[2], point[2]]])

	V1 = np.det(D1)/np.det(D)
	V2 = np.det(D2)/np.det(D)
	V3 = np.det(D3)/np.det(D)
	V  = V1+V2+V3

	return [V1/V, V2/V, V3/V]

##########################################

def sort_split_double(ep1, ep2):
	if abs(ep1[0] - ep2[0]) > 1e-6:
		if ep1[0] - ep2[0] > 0:
			return ep1, ep2 
		elif ep1[0] - ep2[0] < 0:
			return ep2, ep1
	elif abs(ep1[1] - ep2[1]) > 1e-6:
		if ep1[1] - ep2[1] > 0:
			return ep1, ep2 
		elif ep1[1] - ep2[1] < 0:
			return ep1, ep2 

##########################################
		
# def custom_sort(arrays):
    # Define a custom sorting key function
    # def sort_key(array):
        # return array[0], array[1]

    # Sort the arrays based on the custom key function
    # sorted_arrays = sorted(arrays, key=sort_key)

    # return sorted_arrays

def custom_sort(arrays):
	# Convert the list of arrays into a NumPy array
	arrays = np.array(arrays)

	# Sort the arrays based on the first element
	sorted_indices = np.argsort(arrays[:, 0])
	sorted_arrays = arrays[sorted_indices]

	# Check for elements with first elements within 1e-6 of each other
	first_elements = sorted_arrays[:, 0]
	diff = np.abs(np.diff(first_elements))
	split_indices = np.where(diff > 1e-6)[0] + 1

	# Sort the subarrays based on the second element
	final_sorted_arrays = np.split(sorted_arrays, split_indices)
	final_sorted_arrays = [np.array(sorted(subarr, key=lambda x: x[1])) for subarr in final_sorted_arrays]

	# Concatenate the sorted subarrays
	sorted_arrays = np.concatenate(final_sorted_arrays)

	return sorted_arrays.tolist()

##########################################

def write_out(f_single, f_double, f_triple, addr_binodal_pkl, stable_phi):
	print(f"Processing {addr_binodal_pkl}...")
	if os.path.isfile(f_single):
		f_single = open(args.db_single, 'a')
	else:
		f_single = open(args.db_single, 'w')
		if args.ol:
			f_single.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2\n")
		elif args.ww:
			f_single.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3 | w1 | w2 | w3\n")
		else:
			f_single.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3\n")
	
	if os.path.isfile(f_double):
		f_double = open(args.db_double, 'a')
	else:
		f_double = open(args.db_double, 'w')
		if args.ol:
			f_double.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2\n")
		elif args.ww:
			f_double.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3 | w1 | w2 | w3\n")
		else:
			f_double.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3\n")

	if os.path.isfile(f_triple):
		f_triple = open(args.db_triple, 'a')
	else:
		f_triple = open(args.db_triple, 'w')
		if args.ol:
			f_triple.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2\n")
		elif args.ww:
			f_triple.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3 | w1 | w2 | w3\n")
		else:
			f_triple.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3\n")

	g = open(addr_binodal_pkl, 'rb')
	BINODALS = pickle.load(g)
	g.close()

	# find the number of binodals
	num_double = 0
	num_triple = 0
	for idx, f in enumerate(BINODALS["hull_info"]["binodal"]):
		if BINODALS["hull_info"]["function"][idx] is None:
			continue
		elif f[-1] == "two_phase":
			num_double += 1
		elif f[-1] == "three_phase":
			num_triple += 1

	for hull_info in BINODALS["hull_info"]["binodal"]:
		hull_info[0][0], keep = ternary.remove_close_rows(hull_info[0][0], 1e-6)
		hull_info[0][1]       = hull_info[0][1][keep]

	stable, metastable, indices, splits = check_points_for_stability(BINODALS, stable_phi)

	uindices = np.array(np.unique(indices), dtype=int)

	# get the parameters
	chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))
	chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))
	chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))
	vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))
	vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))
	vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", addr_binodal_pkl)).group(1))

	for i in range(len(uindices)):
		mask = (indices == uindices[i])
	
	
	# write the stable points 	
	if args.ol:
		for ri in range(stable.shape[0]):
			f_single.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {stable[ri][0]} | {stable[ri][1]} | 1 | 0 | 0\n")
	elif args.ww:
		for ri in range(stable.shape[0]):
			f_single.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {stable[ri][0]} | {stable[ri][1]} | 1 | 0 | 0 | {stable[ri][0]} | {stable[ri][1]} | 0 | 0 | 0 | 0 | 1 | 0 | 0\n")
	else:
		for ri in range(stable.shape[0]):
			f_single.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {stable[ri][0]} | {stable[ri][1]} | 1 | 0 | 0 | {stable[ri][0]} | {stable[ri][1]} | 0 | 0 | 0 | 0\n")

	
	# collect points per binodal 
	hcount  = -1 
	for hdx, hull_info in enumerate(BINODALS["hull_info"]["binodal"]):
		if hull_info[-1]=="two_phase": 
			hcount      += 1
			mask = (hdx == indices)
			h_metastable = metastable[mask]
			h_idx        = indices   [mask]
			h_splits     = splits    [mask]
			if args.ol:
				for ri in range(h_metastable.shape[0]):
					f_double.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,1]} | 0 | 1 | 0\n")
			elif args.ww:
				for ri in range(h_metastable.shape[0]):
					wts = double_split(h_splits[ri,2*hcount], h_splits[ri, 2*hcount+1], h_metastable[ri])
					f_double.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,1]} | 0 | 1 | 0 | {h_splits[ri,2*hcount,0]} | {h_splits[ri,2*hcount,1]} | {h_splits[ri,2*hcount+1,0]} | {h_splits[ri,2*hcount+1,1]} | 0 | 0 | {wts[0]} | {wts[1]} | 0\n")
			else:
				for ri in range(h_metastable.shape[0]):
					f_double.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,1]} | 0 | 1 | 0 | {h_splits[ri,2*hcount,0]} | {h_splits[ri,2*hcount,1]} | {h_splits[ri,2*hcount+1,0]} | {h_splits[ri,2*hcount+1,1]} | 0 | 0\n")
		
		if hull_info[-1] == "three_phase":
			triangle = custom_sort(hull_info[0])
			# get points 
			mask = (hdx == indices)
			h_metastable = metastable[mask]
			if args.ol:
				for ri in range(h_metastable.shape[0]):
					f_triple.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,1]} | 0 | 0 | 1\n")
			elif args.ww:
				for ri in range(h_metastable.shape[0]):
					wts = triple_split(h_splits[ri,2*hcount], h_splits[ri, 2*hcount+1], h_splits[ri, 2*hcount+2], h_metastable[ri])
					f_triple.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,1]} | 0 | 0 | 1 | {triangle[0][0]} | {triangle[0][1]} | {triangle[1][0]} | {triangle[1][1]} | {triangle[2][0]} | {triangle[2][1]} | {wts[0]} | {wts[1]} | {wts[2]}\n")
			else:
				for ri in range(h_metastable.shape[0]):
					f_triple.write(f"{vs} | {vc} | {vp} | {chisc} | {chips} | {chipc} | {h_metastable[ri,0]} | {h_metastable[ri,0]} | 0 | 0 | 1 | {triangle[0][0]} | {triangle[0][1]} | {triangle[1][0]} | {triangle[1][1]} | {triangle[2][0]} | {triangle[2][1]} \n")
		
	f_single.close()
	f_double.close()
	f_triple.close()

	return 

##########################################

def perpendicular_distance_to_line_segments(point, EP1, EP2):

	# Calculate vector differences between endpoints
	delta = EP2 - EP1

	# Calculate squared magnitudes of each segment
	delta_mag_squared = np.sum(delta**2, axis=1)

	# Avoid division by zero
	mask  = delta_mag_squared != 0.0
	EP1   = EP1[mask]
	EP2   = EP2[mask]
	delta = delta[mask]
	delta_mag_squared = delta_mag_squared[mask]

	# Calculate parametric values (dot products)
	t = np.sum((point - EP1) * delta, axis=1) / delta_mag_squared

	# Compute the closest points on the line segments to the point
	# points1 = EP1[t<=0]
	# points2 = EP2[t>=1]
	neither = np.logical_and(~(t<0), ~(t>1))
	points3 = EP1[neither] + t[neither][:,np.newaxis] * delta[neither] 
	# print(f"1: {len(points1)}, 2: {len(points2)}, 3: {len(points3)}")
	closest_points = points3

	# print(f"len of closest_points = {len(closest_points)}...")

	# Calculate the perpendicular distance
	if len(closest_points) == 0:
		return None, None, None, None 
	else:
		distances = np.linalg.norm(point - closest_points, axis=1)
		return distances[np.argmin(distances)], closest_points[np.argmin(distances)], EP1[np.argmin(distances)], EP2[np.argmin(distances)]

##########################################

def check_points_for_stability(BINODALS, points):
	b_index    = np.zeros(points.shape[0]) - 1
	to_keep    = np.zeros(points.shape[0])
	ntwo_phase = 0

	for hidx, hull in enumerate(BINODALS["hull_info"]["binodal"]):
		if hull[-1] == "two_phase":
			ntwo_phase += 1

	splits     = np.zeros((points.shape[0], 2*ntwo_phase, 2))
	hcount     = -1

	print(f"ntwo_phase = {ntwo_phase}.")

	for hidx, hull in enumerate(BINODALS["hull_info"]["binodal"]):
		if hull[-1] == "two_phase":
			hcount += 1
			h1 = hull[0][0][:,0:2]
			h2 = hull[0][1][:,0:2]
			hh = np.vstack((np.flip(h1[:,0:2], axis=0), h2[:,0:2], h1[-1,0:2]))
			for pidx, pt in enumerate(points):
				# print(f"@ {pidx}/{len(points)}...", flush=True)
				min_distance, closest_point, ep1, ep2 = perpendicular_distance_to_line_segments(pt, h1, h2)
				if min_distance is None:
					to_keep[pidx] += 1
				else:
					eps = custom_sort([ep1, ep2])
					# print(f"eps = {eps}")
					Lh   = LineString(hh)
					dir = (pt[0:2] - closest_point[0:2])/np.linalg.norm(pt[0:2] - closest_point[0:2])
					seg  = np.linspace(pt[0:2]+0.001*dir, closest_point[0:2]-0.001*dir, 1000)
					Lseg = LineString(seg)
					intersection = Lh.intersection(Lseg)

					splits[pidx,2*hcount]   = eps[0]
					splits[pidx,2*hcount+1] = eps[1]
					if intersection.is_empty:
						# this implies that the point is possibly INSIDE a binodal
						if b_index[pidx] == -1:
							b_index[pidx] = hidx
						else:
							print(f"There is an issue. A point is inside two binodals. Exiting...", flush=True)
							exit()
					else:
						to_keep[pidx] += 1
	
	# this is for outside binodal
	net_mask   = (to_keep == ntwo_phase)

	for hidx, hull_path in enumerate(BINODALS["hull_info"]["binodal"]):
		if hull_path[-1] == "three_phase":
			# check if inside binodal
			hpath    = Path(hull_path[0][:,0:2])
			mask     = hpath.contains_points(points[:,0:2], radius=0.001)
			net_mask = np.logical_or(net_mask, ~mask)
			b_index[mask] = hidx
	
	stable     = points [net_mask]
	metastable = points [~net_mask]
	b_index    = b_index[~net_mask]
	splits     = splits [~net_mask]

	# print(f"metastable points are \n{metastable}")
	'''
	for hidx, hull in enumerate(BINODALS["hull_info"]["binodal"]):
		if hull[-1] == "two_phase":
			hcount += 1
			h1 = hull[0][0][:,0:2]
			h2 = hull[0][1][:,0:2]
			hh = np.vstack((np.flip(h1[:,0:2], axis=0), h2[:,0:2], h1[-1,0:2]))
			HULL = ConvexHull(hh)
			HULL = Path(hh[HULL.vertices])
			ax.plot(hh[:,0], 1-hh[:,0]-hh[:,1], hh[:,1], c="gold", lw=0.25)
			for pidx, pt in enumerate(metastable[-1:]):
				print(f"@ {pidx}/{len(metastable)}...", flush=True)
				min_distance, closest_point, ep1, ep2 = perpendicular_distance_to_line_segments(pt, h1, h2)
				eps = custom_sort([ep1, ep2])
				print(f"eps = {eps}")
				Lh   = LineString(hh[:,0:2])
				print(f"is inside: {HULL.contains_point(closest_point[0:2])}")
				print(f"is inside: {HULL.contains_point(pt[0:2])}")
				dir = (pt[0:2] - closest_point[0:2])/np.linalg.norm(pt[0:2] - closest_point[0:2])
				seg  = np.linspace(pt[0:2]+0.001*dir, closest_point[0:2]-0.001*dir, 1000)
				Lseg = LineString(seg)
				# Lseg = LineString(np.vstack((pt[0:2],closest_point[0:2])))
				intersection = Lh.intersection(Lseg)
				print(f"intersection = {intersection}", flush=True)
				intersection = Lseg.intersection(Lh)
				print(f"intersection = {intersection}", flush=True)
				ax.plot(seg[:,0], 1-seg[:,0]-seg[:,1], seg[:,1], c="red", lw=0.25)

				print(f"pt = {pt}: closest point = {closest_point[0:2]}, min dist = {min_distance}.", flush=True)
	'''

	return stable, metastable, b_index, splits

##########################################

def create_uniform_mesh(num_points):
	"""
	Create a uniform mesh in a ternary diagram.

	Parameters:
	- num_points: Number of points along each edge of the triangle

	Returns:
	- mesh_points: Array of mesh points with coordinates (x, y, z)
	"""
	# Initialize lists to store x, y, z coordinates
	x_coords = []
	y_coords = []
	z_coords = []

	# Generate coordinates for the mesh points
	for i in range(num_points + 1):
		for j in range(num_points + 1 - i):
			x = i / num_points
			y = j / num_points
			z = 1 - x - y
			if 0.001 <= x <= 0.999 and 0.001 <= y <= 0.999 and 0.001 <= z <= 0.999:
				x_coords.append(x)
				y_coords.append(y)
				z_coords.append(z)

	# Convert lists to numpy arrays
	x_coords = np.array(x_coords)
	y_coords = np.array(y_coords)
	z_coords = np.array(z_coords)

	# Stack x, y, z coordinates to create mesh points array
	mesh_points = np.column_stack((x_coords, y_coords, z_coords))

	return mesh_points

##########################################

if __name__=="__main__":

	start = time.time()

	#####################################################
	# get a mesh of volume fractions
	mesh               = create_uniform_mesh(args.meshd)
	stable_phi         = mesh[:,0:2]
	
	f_single = args.db_single
	f_double = args.db_double
	f_triple = args.db_triple
	
	for directory in args.dirs:
		l_files = os.listdir(directory)
		for fyle in l_files:
			write_out(f_single, f_double, f_triple, directory+"/"+fyle, stable_phi)

	if args.combine:
		file_paths = [f_single, f_double, f_triple]
		dfs = []
		for file_path in file_paths:
			df = pd.read_csv(file_path, sep='|', header=0) 
			dfs.append(df)

			# Concatenate the dataframes vertically (along rows)
			combined_df = pd.concat(dfs, ignore_index=True)

			# Write the combined dataframe to a new CSV file
			combined_df.to_csv(args.db_comb, sep='|', index=False)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

