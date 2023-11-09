#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import mpltern

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--no-rtw', dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--stable', dest='ons', action='store_true', default=False, help="Only keep fractions that are stable. (default: False).")
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    if args.nrtw:
        return f"beep.\n"
    else:
        return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
#########################################

def generate_points_on_circle(x, r, M):

	# Calculate angles for equally spaced points around the circle
	angles = np.linspace(0, 2 * np.pi, M, endpoint=False)

	# Calculate the coordinates of points on the circle
	x_coordinates = x[0] + r * np.cos(angles)
	y_coordinates = x[1] + r * np.sin(angles)

	points = np.column_stack((x_coordinates, y_coordinates))

	return points

#########################################

def generate_mesh_within_circle(x, rmin, rmax, rdensity, M):

	r_ = np.linspace(rmin, rmax, rdensity)
	for idx,r in enumerate(r_):
		if idx == 0:
			mesh = generate_points_on_circle(x, r, M)
		else:
			mesh = np.vstack((mesh, generate_points_on_circle(x, r, M)))

	return mesh

#########################################

def perform_sweep (chi_ps, chi_pc, chi_sc, crit_point, center):

	# print (f"pid = {os.getpid()}.", flush=True)
	# generate points
	mesh = generate_mesh_within_circle(center, 0.01, 0.1, 100, 50)

	phi_a = mesh[:,0]
	phi_b = mesh[:,1]

	# only keep stuff which is outside the spinodal
	if args.ons:
		to_keep = ternary.stab_crit (phi_a, phi_b, chi_ps, chi_pc, chi_sc) >= 0
		phi_b   = phi_b [to_keep]
		phi_a   = phi_a [to_keep]

	# now start splitting up phi_a, phi_b
	phis    = np.vstack((phi_a, phi_b)).T

	# define the axis along which to make the split
	central_axis   = (crit_point-center) / np.linalg.norm (crit_point-center)
	
	# find those ABOVE axis
	direction      = (phis - center) / np.linalg.norm(phis-center, axis=1)[:, np.newaxis]
	clock          = np.cross (central_axis, direction)
	phi_upper      = phis[clock > 0]
	phi_lower      = phis[clock < 0]

	# now split phi_a, phi_b such that they are on either side of the critical line
	chem_pot_a_upper = mu_a (phi_upper[:,0], phi_upper[:,1])
	chem_pot_b_upper = mu_b (phi_upper[:,0], phi_upper[:,1])
	chem_pot_c_upper = mu_c (phi_upper[:,0], phi_upper[:,1])

	masks = np.isinf(chem_pot_a_upper) | np.isnan(chem_pot_a_upper) | np.isinf(chem_pot_b_upper) |\
	np.isnan(chem_pot_b_upper) | np.isinf (chem_pot_c_upper) | np.isnan (chem_pot_c_upper)

	chem_pot_a_upper = chem_pot_a_upper [~masks]
	chem_pot_b_upper = chem_pot_b_upper [~masks]
	chem_pot_c_upper = chem_pot_c_upper [~masks]

	# get only the relevant phi_uppers
	phi_upper = phi_upper[~masks]

	# calculate the chemical potentials of the lower half
	chem_pot_a_lower = mu_a (phi_lower[:,0], phi_lower[:,1])
	chem_pot_b_lower = mu_b (phi_lower[:,0], phi_lower[:,1])
	chem_pot_c_lower = mu_c (phi_lower[:,0], phi_lower[:,1])

	# create a mask for all the good lower chemical potentials
	masks = np.isinf(chem_pot_a_lower) | np.isnan(chem_pot_a_lower) | np.isinf(chem_pot_b_lower) |\
	np.isnan(chem_pot_b_lower) | np.isinf (chem_pot_c_lower) | np.isnan (chem_pot_c_lower)

	# get only the relevant phi_lowers
	chem_pot_a_lower = chem_pot_a_lower [~masks]
	chem_pot_b_lower = chem_pot_b_lower [~masks]
	chem_pot_c_lower = chem_pot_c_lower [~masks]
	phi_lower = phi_lower[~masks]

	# get an array of these values
	mu_upper   = np.array ([chem_pot_a_upper, chem_pot_b_upper, chem_pot_c_upper]).T
	mu_lower   = np.array ([chem_pot_a_lower, chem_pot_b_lower, chem_pot_c_lower]).T

	# calculate all the distances
	distances       = np.linalg.norm (mu_upper[:, np.newaxis] - mu_lower, axis=2)
	print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
	closest_indices = np.argmin(distances, axis=1)

	# partition, partition... 
	phi_lower       = phi_lower[closest_indices]
	mu_lower        = mu_lower [closest_indices]
	min_distances   = distances[np.arange(len(mu_upper)), closest_indices]

	return [phi_upper, phi_lower, min_distances]

#########################################

def find_binodals(phi_upper, phi_lower, central_axis, center):
	print(f"center = {center}")
	print(f"central_axis = {central_axis}")
	sol_bin_up   = np.empty((0,3))
	sol_bin_down = np.empty((0,3))

	print(f"phi_upper = {phi_upper.shape}")
	for idx, pu in enumerate(phi_upper):
		print(f"idx = {idx}... {idx}/{phi_upper.shape[0]}.")
		def mu_equations(phi):
			eq1 = mu_a(pu[0], phi[0]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(pu[0], phi[0]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(pu[0], phi[0]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for jdx, pl in enumerate(phi_lower):
			root = fsolve(mu_equations, [pu[1], pl[0], pl[1]])
			if (np.abs(np.array(mu_equations(root))>1e-8)).any():
				continue
			else:
				p1 = np.array([pu[0], 1-pu[0]-root[0], root[0]])
				p2 = np.array([root[1], 1-root[1]-root[2], root[2]])

				if np.linalg.norm(p1-p2) < 1e-3:
					continue

				elif np.isnan(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ons and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up,   p1))
						sol_bin_down = np.vstack((sol_bin_down, p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up,   p2))
						sol_bin_down = np.vstack((sol_bin_down, p1))
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)
					break

	for idx, pu in enumerate(phi_upper):
		print(f"idx = {idx}... {idx}/{phi_upper.shape[0]}.")
		def mu_equations(phi):
			eq1 = mu_a(phi[0], pu[1]) - mu_a(phi[1], phi[2])
			eq2 = mu_b(phi[0], pu[1]) - mu_b(phi[1], phi[2])
			eq3 = mu_c(phi[0], pu[1]) - mu_c(phi[1], phi[2])
			return [eq1, eq2, eq3]

		for jdx, pl in enumerate(phi_lower):
			root = fsolve(mu_equations, [pu[1], pl[0], pl[1]])
			if (np.abs(np.array(mu_equations(root))>1e-8)).any():
				continue
			else:
				p1 = np.array([root[0], 1-pu[1]-root[0], pu[1]])
				p2 = np.array([root[1], 1-root[1]-root[2], root[2]])

				if np.linalg.norm(p1-p2) < 1e-3:
					continue

				elif np.isnan(ternary.stab_crit(p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif np.isinf(ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)):
					continue

				elif args.ons and (ternary.stab_crit (p1[0], p1[1], vs, vc, vp, chi_ps, chi_pc, chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], vs, vc, vp, chi_ps, chi_pc, chi_sc) < 0):
					continue

				else:
					if np.sign(np.cross (central_axis, p1[0:2]-center)) == np.sign(np.cross (central_axis, p2[0:2]-center)):
						continue
					elif np.cross (central_axis, p1[0:2]-center)>=0:
						sol_bin_up   = np.vstack((sol_bin_up,   p1))
						sol_bin_down = np.vstack((sol_bin_down, p2))
					else:
						sol_bin_up   = np.vstack((sol_bin_up,   p2))
						sol_bin_down = np.vstack((sol_bin_down, p1))
					print ("HIT!", flush=True, end=' ')
					print (f"p1 = {p1}, p2 = {p2}!", flush=True)
					break

	
	return sol_bin_up, sol_bin_down

#########################################

if __name__=="__main__":

	start = time.time()

	#####################################################
	fig = plt.figure(figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	########################################################
	# set up the inputs
	# set up the chi values
	chi_sc = args.chi_sc
	chi_ps = args.chi_ps
	chi_pc = args.chi_pc

	# set up the volume fractions
	vs = args.vs
	vp = args.vp
	vc = args.vc

	############################################################################################################################
	# solution in terms of phi_s
	discriminant_s = lambda phi_s: -4*vc*vp*(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*\
	vs*chi_pc*(chi_ps+chi_sc))*(phi_s*vs+(-1+phi_s)*vc*(-1+2*phi_s*vs*chi_sc))+(vp - 2*phi_s*vp *vs *chi_ps + \
	vc*(-1+2*phi_s*vs*chi_sc+(-1+phi_s)*vp*(2*chi_pc+phi_s*vs*chi_pc**2 +phi_s*vs*(chi_ps-chi_sc)**2 \
	- 2*phi_s*vs*chi_pc*(chi_ps+chi_sc))))**2

	denom_s    = lambda phi_s: 1/(-2*vc*vp*(2*chi_pc+phi_s*vs*chi_pc**2+phi_s*vs*(chi_ps-chi_sc)**2 - 2*phi_s*vs*chi_pc*(chi_ps+chi_sc)))
	prefac_s   = lambda phi_s: vp - 2*phi_s*vp*vs*chi_ps+vc * (-1+2*phi_s*vs*chi_sc + (-1+phi_s) * vp *\
	(2*chi_pc + phi_s*vs*chi_pc**2 + phi_s * vs * (chi_ps - chi_sc) **2 - 2 * phi_s * vs * chi_pc *(chi_ps + chi_sc)))

	root_up_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) + np.sqrt(discriminant_s(phi_s)))
	root_lo_s  = lambda phi_s: denom_s(phi_s)*(prefac_s(phi_s) - np.sqrt(discriminant_s(phi_s)))

	############################################################################################################################
	# solution in terms of phi_p
	discriminant_p = lambda phi_p: -4*vc*vs*(phi_p*vp+(-1+phi_p)*vc*(-1+2*phi_p*vp*chi_pc))*(2*chi_sc+phi_p*vp*\
	(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc)))+(vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*\
	vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2\
	-2*chi_pc*(1+vs*(chi_ps+chi_sc))))))**2

	denom_p        = lambda phi_p: 1/(-2*vc*vs*(2*chi_sc+phi_p*vp*(chi_pc**2+(chi_ps-chi_sc)**2-2*chi_pc*(chi_ps+chi_sc))))
	prefac_p       = lambda phi_p: vs-2*phi_p*vp*vs*chi_ps+vc*(-1-2*vs*chi_sc+phi_p**2*vp*vs*(chi_pc**2+(chi_ps-chi_sc)**2-\
	2*chi_pc*(chi_ps+chi_sc))+phi_p*(2*vs*chi_sc-vp*(vs*chi_pc**2+vs*(chi_ps-chi_sc)**2-2*chi_pc*(1+vs*(chi_ps+chi_sc)))))

	root_up_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)+np.sqrt(discriminant_p(phi_p)))
	root_lo_p      = lambda phi_p: denom_p(phi_p)*(prefac_p(phi_p)-np.sqrt(discriminant_p(phi_p)))

	#####################################################################

	#####################################################################

	# get all the crit points
	# roots_up, roots_down = ternary.find_crit_point(vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s)

	# put them all together
	# crits      = np.vstack ((roots_up, roots_down))

	# get rid ofthe redundant ones
	# threshold  = 1e-6
	# crits      = ternary.remove_close_rows (crits, threshold)
	crits=np.array([[0.37037037, 0.37037037],[0.37037037,0.25925926],[0.25925926, 0.37037037]])
	print(f"cleaned_crits = \n{crits}")

	tern_b  = True
	edge_b  = True
	crits_b = True
	p_s_space = np.arange(0.001, 1-0.001, 0.001)
	p_s = np.repeat(p_s_space, len(p_s_space))

	p_p = np.zeros(p_s.shape)
	for i in range (len(p_s_space)):
		p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

	vals = ternary.stab_crit (p_s, p_p, vs, vc, vp, chi_ps, chi_pc, chi_sc)

	to_keep = ~np.isnan(vals)

	vals = vals [to_keep]
	p_s  = p_s  [to_keep]
	p_p  = p_p  [to_keep]

	if len(vals) == 0:
		print (f"There will be no critical points and no spinodal region.")

	vmax = np.max(vals)
	vmin = np.min(vals)
	norm = colors.SymLogNorm(0.001, vmin=vmin, vmax=vmax) 
	cols = cm.bwr(norm (vals))

	if np.sign (vmax) == np.sign (vmin):
		if np.sign (vmax) >=0:
			vmin = -vmax
			print (f"There is no unstable region.")
		else:
			vmax = -vmin
			print ("There is mostly unstable region.")

	else:
		print ("there exist unstable regions.")
	
	# plot the thing
	ternary.plot(ax, True, True, True, crits, chi_ps, chi_pc, chi_sc, p_s, p_p, cols, root_up_s, root_lo_s)

	############################################################
	# FIND PHI_B GIVEN PHI_A
	mu_a = lambda phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - vs/vp * phi_b - vs/vc * (1-phi_a-phi_b) + vs * (phi_b**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_sc + phi_b * (1-phi_a-phi_b) * (chi_ps + chi_sc - chi_pc) ) 
	mu_b = lambda phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - vp/vs * phi_a - vp/vc * (1-phi_a-phi_b) + vp * (phi_a**2 * chi_ps + (1-phi_a-phi_b)**2 * chi_pc + phi_a * (1-phi_a-phi_b) * (chi_ps + chi_pc - chi_sc) )
	mu_c = lambda phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - vc/vs * phi_a - vc/vp * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )

	# numerically find the binodal points in the blast radius
	print(f"Computing test binodal points...", flush=True)
	results = perform_sweep(chi_ps, chi_pc, chi_sc, crits[0], np.mean(crits, axis=0))
	print(f"results[0].shape = {results[0].shape}, results[1].shape = {results[1].shape}", flush=True)
	print("Found something!", flush=True)
	distances = results[2]
	# pick a threshold distance
	delta = 1
	mask = distances < delta
	distances = distances[mask]
	phi_arm_1 = results[0][mask]
	phi_arm_2 = results[1][mask]

	##################################################################
	center = np.mean(crits, axis=0)
	central_axis = crits[0]-center
	phi_arm_1, phi_arm_2 = find_binodals(phi_arm_1, phi_arm_2, central_axis, center)

	print(f"phi_arm_1.shape={phi_arm_1.shape}")
	print(f"phi_arm_2.shape={phi_arm_2.shape}")
	ax.scatter(phi_arm_1[:,0], 1-phi_arm_1[:,0]-phi_arm_1[:,1], phi_arm_1[:,1], edgecolors='k', c='steelblue')
	ax.scatter(phi_arm_2[:,0], 1-phi_arm_2[:,0]-phi_arm_2[:,1], phi_arm_2[:,1], edgecolors='k', c='coral')
	###################################################################
	
	# plot the tangent and normal
	ternary.add_tang_norm(ax, False, True, True, crits, vs, vc, vp, chi_pc, chi_ps, chi_sc, root_up_s, root_lo_s)

	# add the plot embellishments
	ternary.embelish(ax, True)
	ax.grid()
	fig.savefig("binodal_blastradius_v0", dpi=1200, bbox_inches="tight")

	stop = time.time()
	print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)
