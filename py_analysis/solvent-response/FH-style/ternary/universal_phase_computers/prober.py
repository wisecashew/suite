import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import tangent
import mpltern
import pickle
import phase
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

import argparse
parser = argparse.ArgumentParser(description='.'         )
parser.add_argument('--chisc',               metavar='chi_sc',  dest='chi_sc',        type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',               metavar='chi_ps',  dest='chi_ps',        type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',               metavar='chi_pc',  dest='chi_pc',        type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',                   metavar='vs',      dest='vs',            type=float, action='store', help='specific volume of solvent.'  )
parser.add_argument('-vc',                   metavar='vc',      dest='vc',            type=float, action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',                   metavar='vp',      dest='vp',            type=float, action='store', help='specific volume of polymer.'  )
parser.add_argument('--database',            metavar='DB',      dest='db',            type=str,   action='store', help='name of database for this system.')
parser.add_argument('--crit-pkl',            metavar='critpkl', dest='critpkl',       type=str,   action='store', help='location of serialized critical point (default: None).',                                default=None)
parser.add_argument('--binodal-pkl',         metavar='bpkl',    dest='bpkl',          type=str,   action='store', help='enter name of file with all the information about the binodals (default: None).',       default=None)
parser.add_argument('--plot-edges',     dest='pe',         action='store_true',  help='plot the edges of the spinodal.',     default=False)
parser.add_argument('--plot-crits',     dest='pc',         action='store_true',  help='plot the critical points.'      ,     default=False)
parser.add_argument('--plot-binodals',  dest='pb',         action='store_true',  help='plot the binodal points.'       ,     default=False)
parser.add_argument('--no-rtw',         dest='nrtw',       action='store_true',  help="Dont print out the runtime warning.", default=False)
parser.add_argument('--img',            dest='img',        action='store',       type=str, default="None", help="name of image to be created.")
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	if args.nrtw:
		return f""
	else:
		return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format

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
	return V1/V, V2/V

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

	return V1/V, V2/V, V3/V 

##########################################

if __name__=="__main__":

	start = time.time()

	#####################################################
	fig = plt.figure(num=0, figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	########################################################
	# set up the inputs
	# set up the chi values
	inputs = dict()
	inputs["chi_sc"] = args.chi_sc
	inputs["chi_ps"] = args.chi_ps
	inputs["chi_pc"] = args.chi_pc
	inputs["vs"] = args.vs
	inputs["vp"] = args.vp
	inputs["vc"] = args.vc
	tern_b  = True
	edges_b = args.pe
	crits_b = args.pc

	# set up objects and crit points
	print(f"Setting up objects...", flush=True, end=' ')
	P = phase.Phase(inputs)
	print("done!", flush=True)

	if args.critpkl == None:	
		P.spinodal.obtain_crits()
	else:
		try:
			f = open(args.critpkl, 'rb')
			crits = pickle.load(f)
			P.spinodal.crits = crits
			f.close()
		except:
			P.spinodal.obtain_crits()
			f = open(args.critpkl, 'wb')
			pickle.dump(P.spinodal.crits, f)
			f.close()
	
	P.crits = P.spinodal.crits
	print(f"crits = {P.crits}", flush=True)
	print(f"done!", flush=True)

	print(f"Plotting the ternary diagram...", flush=True,end=' ')
	P.spinodal.stability_plots(ax, tern_b, edges_b, crits_b)
	print(f"done!", flush=True)

	# P.tangent_tracing(ax)
	# print("Plotted out the tangent trace!", flush=True)
	
	# get a mesh of volume fractions
	phi_s, phi_p = np.meshgrid(np.linspace(0.01,0.99,50), np.linspace(0.01,0.99,50))
	mask     = (phi_s + phi_p < 1)
	phi_s    = phi_s[mask] 
	phi_p    = phi_p[mask]
	stable_phi = np.array([phi_s, phi_p]).T

	f = open(args.bpkl, 'rb')
	BINODALS = pickle.load(f)

	for hidx,hull in enumerate(BINODALS["hull_info"]["function"]):
		print(f'arm 1 = {BINODALS["hull_info"]["binodal"][hidx][0][0][-10:]}')
		print(f'arm 2 = {BINODALS["hull_info"]["binodal"][hidx][0][1][-10:]}')
		mask       = hull.contains_points(stable_phi) 
		stable_phi = stable_phi[~mask]
	
	ax.scatter(stable_phi[:,0], 1-stable_phi[:,0]-stable_phi[:,1], stable_phi[:,1], c='lavender')

	f = open(args.db, 'w')
	f.write(f"vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3\n")
	for sphi in stable_phi:
		f.write(f"{P.vs} | {P.vc} | {P.vp} | {P.chi_sc} | {P.chi_ps} | {P.chi_pc} | {sphi[0]} | {sphi[1]} | {sphi[0]} | {sphi[1]} | 0 | 0 | 0 | 0\n")


	for hull_info in BINODALS["hull_info"]["binodal"]:
		if hull_info[-1]=="two_phase":
			# print(f"arm 1 = {hull_info[0][0][-10:]}")
			# print(f"arm 2 = {hull_info[0][1][-10:]}")
			if len(hull_info) == 2:
				ax.scatter(hull_info[0][0][:,0], 1-hull_info[0][0][:,0]-hull_info[0][0][:,1], hull_info[0][0][:,1], c='black', s=0.5)
				ax.scatter(hull_info[0][1][:,0], 1-hull_info[0][1][:,0]-hull_info[0][1][:,1], hull_info[0][1][:,1], c='white', s=0.5)
				for i in range(len(hull_info[0][0])):
					line = np.linspace(hull_info[0][0][i][0:2], hull_info[0][1][i][0:2],100)
					ax.scatter(line[:,0], 1-line[:,0]-line[:,1], line[:,1], c='pink', s=0.5)
					for p in line:
						f.write(f"{P.vs} | {P.vc} | {P.vp} | {P.chi_sc} | {P.chi_ps} | {P.chi_pc} | {p[0]} | {p[1]} | {line[0][0]} | {line[0][1]} | {line[-1][0]} | {line[-1][1]} | 0 | 0\n")

			elif len(hull_info) == 3:
				ax.scatter(hull_info[0][:,0], 1-hull_info[0][:,0]-hull_info[0][:,1], hull_info[0][:,1], c='black', s=0.5)
				ax.scatter(hull_info[1][:,0], 1-hull_info[1][:,0]-hull_info[1][:,1], hull_info[1][:,1], c='white', s=0.5)
				for i in range(len(hull_info[0])):
					line = np.linspace(hull_info[0][i][0:2], hull_info[1][i][0:2],100)
					ax.scatter(line[:,0], 1-line[:,0]-line[:,1], line[:,1], c='pink', s=0.5)
					for p in line:
						f.write(f"{P.vs} | {P.vc} | {P.vp} | {P.chi_sc} | {P.chi_ps} | {P.chi_pc} | {p[0]} | {p[1]} | {line[0][0]} | {line[0][1]} | {line[-1][0]} | {line[-1][1]} | 0 | 0\n")

		elif hull_info[-1]=="three_phase":
			# get points 
			in_triangle = generate_points_inside_triangle(hull_info[0], 1000)
			ax.scatter(in_triangle[:,0], 1-in_triangle[:,0]-in_triangle[:,1], in_triangle[:,1], c='lawngreen', s=0.5)
			tri_pad = np.vstack((hull_info[0], hull_info[0][0]))
			print(f"triangle = {tri_pad}")
			ax.plot(tri_pad[:,0], 1-tri_pad[:,0]-tri_pad[:,1], tri_pad[:,1], c='darkred', lw=0.5)
			for p in in_triangle:
				f.write(f"{P.vs} | {P.vc} | {P.vp} | {P.chi_sc} | {P.chi_ps} | {P.chi_pc} | {p[0]} | {p[1]} | {tri_pad[0][0]} | {tri_pad[0][1]} | {tri_pad[1][0]} | {tri_pad[1][1]} | {tri_pad[2][0]} | {tri_pad[2][1]}\n")
	f.close()

	# create the image
	print("Making image...", end=' ', flush=True)
	
	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		fig.savefig (img_name, dpi=1200, bbox_inches="tight")
	else:
		fig.savefig (f"internal_points-vs_{P.vs}-vc_{P.vc}-vp_{P.vp}-chisc_{P.chi_sc}-chips_{P.chi_ps}-chipc_{P.chi_pc}.png", dpi=1200)
	
	print(f"done!", flush=True)

	stop = time.time()
	print(f"Time for computation is {stop-start} seconds.", flush=True)

