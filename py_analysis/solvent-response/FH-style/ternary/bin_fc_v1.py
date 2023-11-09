#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.optimize import fsolve
import argparse
import time
import warnings
import linecache
import ternary
import tangent
import mpltern

import argparse
parser = argparse.ArgumentParser(description="Create a skeleton solution for the binodal. This is a memory-intensive computation.")
parser.add_argument('--chisc',  metavar='chi_sc', dest='chi_sc',  type=float, action='store', help='enter S-C exchange parameter.')
parser.add_argument('--chips',  metavar='chi_ps', dest='chi_ps',  type=float, action='store', help='enter P-S exchange parameter.')
parser.add_argument('--chipc',  metavar='chi_pc', dest='chi_pc',  type=float, action='store', help='enter P-C exchange parameter.')
parser.add_argument('-vs',      metavar='vs',     dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',     dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',     dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--iter',   metavar='it',     dest='it',      type=int,     action='store', help='Specify the number of secondary searches you want (default: 50).', default=50)
parser.add_argument('--no-rtw', dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--stable', dest='ons', action='store_true', default=False, help="Only keep fractions that are stable. (default: False).")
parser.add_argument('--img',    dest='img', action='store',      type=str, default="None", help="name of image to be created. (default: blastradius).")
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

	row_sums = np.sum(mesh, axis=1)
	mask     = (row_sums <= 1)
	mesh     = mesh[mask]

	return mesh

#########################################

#########################################

#########################################

class Binodal:

	def __init__(self, inputs, crits=None):
		self.chi_sc = inputs["chi_sc"]
		self.chi_ps = inputs["chi_ps"]
		self.chi_pc = inputs["chi_pc"]
		self.vs     = inputs["vs"]
		self.vc     = inputs["vc"]
		self.vp     = inputs["vp"]
		self.crits  = crits 
		return
	
	# solution in terms of phi_s
	discriminant_s = lambda self,phi_s: -4*self.vc*self.vp*(2*self.chi_pc + phi_s*self.vs*self.chi_pc**2 + phi_s*self.vs*(self.chi_ps-self.chi_sc)**2 - 2*phi_s*\
	self.vs*self.chi_pc*(self.chi_ps+self.chi_sc))*(phi_s*self.vs+(-1+phi_s)*self.vc*(-1+2*phi_s*self.vs*self.chi_sc))+(self.vp - 2*phi_s*self.vp*self.vs*self.chi_ps + \
	self.vc*(-1+2*phi_s*self.vs*self.chi_sc+(-1+phi_s)*self.vp*(2*self.chi_pc+phi_s*self.vs*self.chi_pc**2 +phi_s*self.vs*(self.chi_ps-self.chi_sc)**2 \
	- 2*phi_s*self.vs*self.chi_pc*(self.chi_ps+self.chi_sc))))**2

	denom_s    = lambda self,phi_s: 1/(-2*self.vc*self.vp*(2*self.chi_pc+phi_s*self.vs*self.chi_pc**2+phi_s*self.vs*(self.chi_ps-self.chi_sc)**2 - 2*phi_s*self.vs*self.chi_pc*(self.chi_ps+self.chi_sc)))
	prefac_s   = lambda self,phi_s: self.vp - 2*phi_s*self.vp*self.vs*self.chi_ps+self.vc * (-1+2*phi_s*self.vs*self.chi_sc + (-1+phi_s) * self.vp *\
	(2*self.chi_pc + phi_s*self.vs*self.chi_pc**2 + phi_s * self.vs * (self.chi_ps - self.chi_sc) **2 - 2 * phi_s * self.vs * self.chi_pc *(self.chi_ps + self.chi_sc)))

	root_up_s  = lambda self,phi_s: self.denom_s(phi_s)*(self.prefac_s(phi_s) + np.sqrt(self.discriminant_s(phi_s)))
	root_lo_s  = lambda self,phi_s: self.denom_s(phi_s)*(self.prefac_s(phi_s) - np.sqrt(self.discriminant_s(phi_s)))

	############################################################################################################################
	# solution in terms of phi_p
	discriminant_p = lambda self,phi_p: -4*self.vc*self.vs*(phi_p*self.vp+(-1+phi_p)*self.vc*(-1+2*phi_p*self.vp*self.chi_pc))*(2*self.chi_sc+phi_p*self.vp*\
	(self.chi_pc**2+(self.chi_ps-self.chi_sc)**2-2*self.chi_pc*(self.chi_ps+self.chi_sc)))+(self.vs-2*phi_p*self.vp*self.vs*self.chi_ps+self.vc*(-1-2*self.vs*self.chi_sc+phi_p**2*\
	self.vp*self.vs*(self.chi_pc**2+(self.chi_ps-self.chi_sc)**2-2*self.chi_pc*(self.chi_ps+self.chi_sc))+phi_p*(2*self.vs*self.chi_sc-self.vp*(self.vs*self.chi_pc**2+self.vs*(self.chi_ps-self.chi_sc)**2\
	-2*self.chi_pc*(1+self.vs*(self.chi_ps+self.chi_sc))))))**2

	denom_p        = lambda self,phi_p: 1/(-2*self.vc*self.vs*(2*self.chi_sc+phi_p*self.vp*(self.chi_pc**2+(self.chi_ps-self.chi_sc)**2-2*self.chi_pc*(self.chi_ps+self.chi_sc))))
	prefac_p       = lambda self,phi_p: self.vs-2*phi_p*self.vp*self.vs*self.chi_ps+self.vc*(-1-2*self.vs*self.chi_sc+phi_p**2*self.vp*self.vs*(self.chi_pc**2+(self.chi_ps-self.chi_sc)**2-\
	2*self.chi_pc*(self.chi_ps+self.chi_sc))+phi_p*(2*self.vs*self.chi_sc-self.vp*(self.vs*self.chi_pc**2+self.vs*(self.chi_ps-self.chi_sc)**2-2*self.chi_pc*(1+self.vs*(self.chi_ps+self.chi_sc)))))

	root_up_p      = lambda self,phi_p: self.denom_p(phi_p)*(self.prefac_p(phi_p)+np.sqrt(self.discriminant_p(phi_p)))
	root_lo_p      = lambda self,phi_p: self.denom_p(phi_p)*(self.prefac_p(phi_p)-np.sqrt(self.discriminant_p(phi_p)))

	# calculate chemical potentials
	mu_a = lambda self, phi_a, phi_b: np.log(phi_a)         + 1 - phi_a - self.vs/self.vp * phi_b - self.vs/self.vc * (1-phi_a-phi_b) + self.vs * (phi_b**2 * self.chi_ps + (1-phi_a-phi_b)**2 * self.chi_sc + phi_b * (1-phi_a-phi_b) * (self.chi_ps + self.chi_sc - self.chi_pc) ) 
	mu_b = lambda self, phi_a, phi_b: np.log(phi_b)         + 1 - phi_b - self.vp/self.vs * phi_a - self.vp/self.vc * (1-phi_a-phi_b) + self.vp * (phi_a**2 * self.chi_ps + (1-phi_a-phi_b)**2 * self.chi_pc + phi_a * (1-phi_a-phi_b) * (self.chi_ps + self.chi_pc - self.chi_sc) )
	mu_c = lambda self, phi_a, phi_b: np.log(1-phi_a-phi_b) + 1 - (1-phi_a-phi_b) - self.vc/self.vs * phi_a - self.vc/self.vp * phi_b + self.vc * (phi_a**2 * self.chi_sc + phi_b**2 * self.chi_pc + phi_a * phi_b * (self.chi_sc + self.chi_pc - self.chi_ps) )

	# get all the crit points
	def obtain_crits(self):
		roots_up, roots_down = ternary.find_crit_point(self.vs, self.vc, self.vp, self.chi_sc, self.chi_ps, self.chi_pc, self.root_up_p, self.root_up_s, self.root_lo_p, self.root_lo_s)
		crits      = np.vstack ((roots_up, roots_down))

		# get rid of the redundant ones
		threshold  = 1e-6
		crits      = ternary.remove_close_rows (crits, threshold)
		self.crits = crits
		return
	######################################################
	# plot the stability criterion on the ternary plot

	def stability_plots(self, ax):
		p_s_space = np.arange(0.001, 1-0.001, 0.001)
		p_s = np.repeat(p_s_space, len(p_s_space))

		p_p = np.zeros(p_s.shape)
		for i in range (len(p_s_space)):
			p_p[i*len(p_s_space):(i+1)*len(p_s_space)] = np.linspace (0.001, 1-p_s_space[i], len(p_s_space))

		vals = ternary.stab_crit (p_s, p_p, self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)

		to_keep = ~np.isnan(vals)

		vals = vals [to_keep]
		p_s  = p_s  [to_keep]
		p_p  = p_p  [to_keep]

		if len(vals) == 0:
			print (f"There will be no critical points and no spinodal region.", flush=True)

		vmax = np.max(vals)
		vmin = np.min(vals)
		norm = colors.SymLogNorm(0.001, vmin=vmin, vmax=vmax) 
		cols = cm.bwr(norm (vals))

		if np.sign (vmax) == np.sign (vmin):
			if np.sign (vmax) >=0:
				vmin = -vmax
				print (f"There is no unstable region.", flush=True)
			else:
				vmax = -vmin
				print ("There is mostly unstable region.", flush=True)

		else:
			print ("there exist unstable regions.", flush=True)
		
		# plot the thing
		ternary.plot(ax, True, True, True, self.crits, self.chi_ps, self.chi_pc, self.chi_sc, p_s, p_p, cols, self.root_up_s, self.root_lo_s)
		ternary.embelish(ax, True)
		ax.grid()
		return
	######################################################
	# start performing sweeps
	# first sweep

	def perform_first_sweep (self, center, crit_point, p_gen):

		# generate points
		mesh = generate_mesh_within_circle(center, p_gen["rmin"], p_gen["rmax"], p_gen["r_dens"], p_gen["theta_dens"])
		phi_s = mesh[:,0]
		phi_p = mesh[:,1]

		# decide if you only keep stuff which is outside the spinodal
		if args.ons:
			to_keep = ternary.stab_crit (phi_s, phi_p, self.chi_ps, self.chi_pc, self.chi_sc) >= 0
			phi_p   = phi_p [to_keep]
			phi_s   = phi_s [to_keep]

		# now start splitting up phi_a, phi_b
		phis    = np.vstack((phi_s, phi_p)).T

		# define the axis along which to make the split
		central_axis   = (crit_point-center) / np.linalg.norm (crit_point-center)
		
		# find those ABOVE axis
		direction      = (phis - center) / np.linalg.norm(phis-center, axis=1)[:, np.newaxis]
		clock          = np.cross (central_axis, direction)
		phi_positive   = phis[clock > 0]
		phi_negative   = phis[clock < 0]

		# now split phi_a, phi_b such that they are on either side of the critical line
		chem_pot_a_positive = self.mu_a(phi_positive[:,0], phi_positive[:,1])
		chem_pot_b_positive = self.mu_b(phi_positive[:,0], phi_positive[:,1])
		chem_pot_c_positive = self.mu_c(phi_positive[:,0], phi_positive[:,1])

		masks = np.isinf(chem_pot_a_positive) | np.isnan(chem_pot_a_positive) | np.isinf(chem_pot_b_positive) |\
		np.isnan(chem_pot_b_positive) | np.isinf (chem_pot_c_positive) | np.isnan (chem_pot_c_positive)

		chem_pot_a_positive = chem_pot_a_positive[~masks]
		chem_pot_b_positive = chem_pot_b_positive[~masks]
		chem_pot_c_positive = chem_pot_c_positive[~masks]

		# get only the relevant phi_uppers
		phi_positive = phi_positive[~masks]

		# calculate the chemical potentials of the lower half
		chem_pot_a_negative = self.mu_a(phi_negative[:,0], phi_negative[:,1])
		chem_pot_b_negative = self.mu_b(phi_negative[:,0], phi_negative[:,1])
		chem_pot_c_negative = self.mu_c(phi_negative[:,0], phi_negative[:,1])

		# create a mask for all the good lower chemical potentials
		masks = np.isinf(chem_pot_a_negative) | np.isnan(chem_pot_a_negative) | np.isinf(chem_pot_b_negative) |\
		np.isnan(chem_pot_b_negative) | np.isinf (chem_pot_c_negative) | np.isnan (chem_pot_c_negative)

		# get only the relevant phi_lowers
		chem_pot_a_negative = chem_pot_a_negative[~masks]
		chem_pot_b_negative = chem_pot_b_negative[~masks]
		chem_pot_c_negative = chem_pot_c_negative[~masks]
		phi_negative = phi_negative[~masks]

		# get an array of these values
		mu_positive   = np.array ([chem_pot_a_positive, chem_pot_b_positive, chem_pot_c_positive]).T
		mu_negative   = np.array ([chem_pot_a_negative, chem_pot_b_negative, chem_pot_c_negative]).T

		# calculate all the distances
		distances       = np.linalg.norm (mu_positive[:, np.newaxis] - mu_negative, axis=2)
		print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
		closest_indices = np.argmin(distances, axis=1)

		# partition, partition... 
		# phi_positive       = phi_positive[closest_indices]
		phi_negative       = phi_negative[closest_indices]
		min_distances      = distances[np.arange(len(mu_positive)), closest_indices]

		return phi_positive, phi_negative, min_distances
	######################################################

	def perform_first_binodal_blast (self, phi_upper, phi_lower, center, crit_point):
		central_axis = (crit_point - center)/np.linalg.norm(crit_point - center)
		sol_bin_up   = np.empty((0,2))
		sol_bin_down = np.empty((0,2))
		print(f"First half...", flush=True)
		for idx, pu in enumerate(phi_upper):
			if idx%100 == 0:
				print(f"First half @ idx = {idx}/{phi_upper.shape[0]}", flush=True)
			def mu_equations(phi):
				eq1 = self.mu_a(pu[0], phi[0]) - self.mu_a(phi[1], phi[2])
				eq2 = self.mu_b(pu[0], phi[0]) - self.mu_b(phi[1], phi[2])
				eq3 = self.mu_c(pu[0], phi[0]) - self.mu_c(phi[1], phi[2])
				return [eq1, eq2, eq3]
			
			for jdx, pl in enumerate(phi_lower):
				root = fsolve(mu_equations, [pu[1], pl[0], pl[1]], xtol=1e-14)
				if (np.abs(np.array(mu_equations(root)))>1e-13).any():
					continue
				else:
					p1 = np.array([pu[0], root[0]])
					p2 = np.array([root[1], root[2]])				
					if np.linalg.norm(p1-p2) < 1e-3:
						continue

					elif np.isnan(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif np.isinf(ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif args.ons and (ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0):
						continue

					else:
						if np.sign(np.cross (central_axis, p1-center)) == np.sign(np.cross (central_axis, p2-center)):
							continue
						elif np.cross (central_axis, p1[0:2]-center)>=0:
							sol_bin_up   = np.vstack((sol_bin_up,   p1))
							sol_bin_down = np.vstack((sol_bin_down, p2))
						else:
							sol_bin_up   = np.vstack((sol_bin_up,   p2))
							sol_bin_down = np.vstack((sol_bin_down, p1))
						print (f"HIT! ({idx}/{phi_upper.shape[0]})", flush=True, end=' ')
						print (f"jdx = {jdx}, p1 = {p1}, p2 = {p2}!", flush=True)
						print(f"seed was phi_s = {pu[0]}, guess was phi_p = {pu[1]}, phi = {pl[0], pl[1]}", flush=True)
						# print(f"delta_mu_init = {mu_equations([pu[1], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						if abs(p1[0]-p2[1])>1e-3 or abs(p1[1]-p2[0])>1e-3:
							print(f"There is a problem in the first half of loop in primary...", flush=True)
							print(f"seed was phi_s = {pu[0]}, guess was phi_p = {pu[1]}, phi = {pl[0], pl[1]}", flush=True)
							print(f"delta_mu_init = {mu_equations([pu[1], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}", flush=True)
							
						break
		
		print(f"Second half...", flush=True)
		for idx, pu in enumerate(phi_upper):
			if idx% 100 == 0:
				print(f"@ idx = {idx}/{phi_upper.shape[0]}", flush=True)
			def mu_equations(phi):
				eq1 = self.mu_a(phi[0], pu[1]) - self.mu_a(phi[1], phi[2])
				eq2 = self.mu_b(phi[0], pu[1]) - self.mu_b(phi[1], phi[2])
				eq3 = self.mu_c(phi[0], pu[1]) - self.mu_c(phi[1], phi[2])
				return [eq1, eq2, eq3]

			for jdx, pl in enumerate(phi_lower):
				root = fsolve(mu_equations, [pu[0], pl[0], pl[1]], xtol=1e-14)
				if (np.abs(np.array(mu_equations(root)))>1e-13).any():
					continue
				else:
					p1 = np.array([root[0], pu[1]])
					p2 = np.array([root[1], root[2]])

					if np.linalg.norm(p1-p2) < 1e-3:
						continue

					elif np.isnan(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif np.isinf(ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif args.ons and (ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0):
						continue

					else:
						if np.sign(np.cross (central_axis, p1-center)) == np.sign(np.cross (central_axis, p2-center)):
							continue
						elif np.cross (central_axis, p1-center)>=0:
							sol_bin_up   = np.vstack((sol_bin_up,   p1))
							sol_bin_down = np.vstack((sol_bin_down, p2))
						else:
							sol_bin_up   = np.vstack((sol_bin_up,   p2))
							sol_bin_down = np.vstack((sol_bin_down, p1))
						print (f"HIT! ({idx}/{phi_upper.shape[0]})", flush=True, end=' ')
						print (f"jdx = {jdx}, p1 = {p1}, p2 = {p2}!", flush=True)
						# print(f"delta_mu_init = {mu_equations([pu[0], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						if abs(p1[0]-p2[1])>1e-3 or abs(p1[1]-p2[0])>1e-3:
							print(f"There is a problem in the second half of loop in primary...", flush=True)
							print(f"seed was phi_p = {pu[1]}, guess was phi_s = {pu[0]}, phi = {pl[0], pl[1]}", flush=True)
							print(f"delta_mu_init = {mu_equations([pu[0], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}", flush=True)
							
						break

		return sol_bin_up, sol_bin_down
	######################################################

	def perform_secondary_sweeps(self, positive_arm_tip, negative_arm_tip, center, crit_point, p_gen):
		print("Performing secondary sweeps...", flush=True)
		print(f"p_gen['rmin'] = {p_gen['rmin']}, p_gen['rmax'] = {p_gen['rmax']}...", flush=True)
		central_axis = (crit_point - center)/np.linalg.norm(crit_point-center)

		# generate points
		mesh_positive  = generate_mesh_within_circle(positive_arm_tip, p_gen["rmin"], p_gen["rmax"], p_gen["r_dens"], p_gen["theta_dens"])
		mesh_negative  = generate_mesh_within_circle(negative_arm_tip, p_gen["rmin"], p_gen["rmax"], p_gen["r_dens"], p_gen["theta_dens"])
		phi_s_positive = mesh_positive[:,0]
		phi_p_positive = mesh_positive[:,1]

		phi_s_negative = mesh_negative[:,0]
		phi_p_negative = mesh_negative[:,1]

		# decide if you only keep stuff which is outside the spinodal
		if args.ons:
			to_keep = ternary.stab_crit(phi_s_positive, phi_p_positive, self.chi_ps, self.chi_pc, self.chi_sc) >= 0
			phi_p_positive   = phi_p_positive[to_keep]
			phi_s_positive   = phi_s_positive[to_keep]
			to_keep = ternary.stab_crit(phi_s_negative, phi_p_negative, self.chi_ps, self.chi_pc, self.chi_sc) >= 0
			phi_p_negative   = phi_p_negative[to_keep]
			phi_s_negative   = phi_s_negative[to_keep]

		# now start splitting up phi_a, phi_b
		phi_positive    = np.vstack((phi_s_positive, phi_p_positive)).T
		phi_negative    = np.vstack((phi_s_negative, phi_p_negative)).T

		# now split phi_a, phi_b such that they are on either side of the critical line
		chem_pot_a_positive = self.mu_a (phi_positive[:,0], phi_positive[:,1])
		chem_pot_b_positive = self.mu_b (phi_positive[:,0], phi_positive[:,1])
		chem_pot_c_positive = self.mu_c (phi_positive[:,0], phi_positive[:,1])

		# create a mask for all the good upper chemical potentials
		masks_pos = np.isinf(chem_pot_a_positive) | np.isnan(chem_pot_a_positive) | np.isinf(chem_pot_b_positive) |\
		np.isnan(chem_pot_b_positive) | np.isinf (chem_pot_c_positive) | np.isnan (chem_pot_c_positive)

		# calculate the chemical potentials of the lower half
		chem_pot_a_negative = self.mu_a (phi_negative[:,0], phi_negative[:,1])
		chem_pot_b_negative = self.mu_b (phi_negative[:,0], phi_negative[:,1])
		chem_pot_c_negative = self.mu_c (phi_negative[:,0], phi_negative[:,1])

		# create a mask for all the good lower chemical potentials
		masks_neg = np.isinf(chem_pot_a_negative) | np.isnan(chem_pot_a_negative) | np.isinf(chem_pot_b_negative) |\
		np.isnan(chem_pot_b_negative) | np.isinf (chem_pot_c_negative) | np.isnan (chem_pot_c_negative)

		chem_pot_a_positive = chem_pot_a_positive[~masks_pos]
		chem_pot_b_positive = chem_pot_b_positive[~masks_pos]
		chem_pot_c_positive = chem_pot_c_positive[~masks_pos]

		# get only the relevant phi_uppers
		phi_positive = phi_positive[~masks_pos]

		# get only the relevant phi_lowers
		chem_pot_a_negative = chem_pot_a_negative [~masks_neg]
		chem_pot_b_negative = chem_pot_b_negative [~masks_neg]
		chem_pot_c_negative = chem_pot_c_negative [~masks_neg]

		# get only the relevant phi_lowers
		phi_negative = phi_negative[~masks_neg]

		# get an array of these values
		mu_positive   = np.array ([chem_pot_a_positive, chem_pot_b_positive, chem_pot_c_positive]).T
		mu_negative   = np.array ([chem_pot_a_negative, chem_pot_b_negative, chem_pot_c_negative]).T

		print(f"mu_positive = {mu_positive.shape}", flush=True)
		print(f"mu_negative = {mu_negative.shape}", flush=True)

		# calculate all the distances
		distances       = np.linalg.norm (mu_positive[:, np.newaxis] - mu_negative, axis=2)
		print (f"Memory size of distances is {distances.nbytes/1e+9} gigabytes.", flush=True)
		closest_indices = np.argmin(distances, axis=1)

		# partition, partition... 
		# phi_positive       = phi_positive[closest_indices]
		phi_negative       = phi_negative[closest_indices]
		min_distances      = distances[np.arange(len(mu_positive)), closest_indices]

		return phi_positive, phi_negative, min_distances
	######################################################

	def perform_secondary_binodal_blast(self, phi_positive, phi_negative, center, crit_point):
		central_axis = (crit_point - center)/np.linalg.norm(crit_point - center)
		sol_bin_up   = np.empty((0,2))
		sol_bin_down = np.empty((0,2))

		for idx, pu in enumerate(phi_positive):
			if idx%100==0:
				print(f"idx = {idx}, {idx}/{phi_positive.shape[0]}...", flush=True)
			def mu_equations(phi):
				eq1 = self.mu_a(pu[0], phi[0]) - self.mu_a(phi[1], phi[2])
				eq2 = self.mu_b(pu[0], phi[0]) - self.mu_b(phi[1], phi[2])
				eq3 = self.mu_c(pu[0], phi[0]) - self.mu_c(phi[1], phi[2])
				return [eq1, eq2, eq3]

			for jdx, pl in enumerate(phi_negative):
				root = fsolve(mu_equations, [pu[1], pl[0], pl[1]], xtol=1e-14)
				if (np.abs(np.array(mu_equations(root)))>1e-13).any():
					continue
				else:
					p1 = np.array([pu[0],   root[0]])
					p2 = np.array([root[1], root[2]])

					if np.linalg.norm(p1-p2) < 1e-3:
						continue

					elif np.isnan(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif np.isinf(ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif args.ons and (ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0):
						continue

					else:
						if np.sign(np.cross (central_axis, p1-center)) == np.sign(np.cross (central_axis, p2-center)):
							continue
						elif np.cross (central_axis, p1-center)>=0:
							sol_bin_up   = np.vstack((sol_bin_up,   p1))
							sol_bin_down = np.vstack((sol_bin_down, p2))
						else:
							sol_bin_up   = np.vstack((sol_bin_up,   p2))
							sol_bin_down = np.vstack((sol_bin_down, p1))
						print ("HIT!", flush=True, end=' ')
						print (f"jdx = {jdx}, p1 = {p1}, p2 = {p2}!", flush=True)
						# print(f"delta_mu_init = {mu_equations([pu[1], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						if abs(p1[0]-p2[1])>1e-3 or abs(p1[1]-p2[0])>1e-3:
							print(f"There is a problem in the first half of loop in secondary...", flush=True)
							print(f"seed was phi_s = {pu[0]}, guess was phi_p = {pu[1]}, phi = {pl[0], pl[1]}", flush=True)
							print(f"delta_mu_init = {mu_equations([pu[1], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						break

		for idx, pu in enumerate(phi_positive):
			if idx%100==0:
				print(f"idx = {idx}, {idx}/{phi_positive.shape[0]}...", flush=True)
			def mu_equations(phi):
				eq1 = self.mu_a(phi[0], pu[1]) - self.mu_a(phi[1], phi[2])
				eq2 = self.mu_b(phi[0], pu[1]) - self.mu_b(phi[1], phi[2])
				eq3 = self.mu_c(phi[0], pu[1]) - self.mu_c(phi[1], phi[2])
				return [eq1, eq2, eq3]

			for jdx, pl in enumerate(phi_negative):
				root = fsolve(mu_equations, [pu[0], pl[0], pl[1]], xtol=1e-14)
				if (np.abs(np.array(mu_equations(root)))>1e-13).any():
					continue
				else:
					p1 = np.array([root[0], pu[1]])
					p2 = np.array([root[1], root[2]])

					if np.linalg.norm(p1-p2) < 1e-3:
						continue

					elif np.isnan(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif np.isinf(ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue

					elif args.ons and (ternary.stab_crit (p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0):
						continue

					else:
						if np.sign(np.cross (central_axis, p1-center)) == np.sign(np.cross (central_axis, p2-center)):
							continue
						elif np.cross (central_axis, p1-center)>=0:
							sol_bin_up   = np.vstack((sol_bin_up,   p1))
							sol_bin_down = np.vstack((sol_bin_down, p2))
						else:
							sol_bin_up   = np.vstack((sol_bin_up,   p2))
							sol_bin_down = np.vstack((sol_bin_down, p1))
						print ("HIT!", flush=True, end=' ')
						print (f"jdx = {jdx}, p1 = {p1}, p2 = {p2}!", flush=True)
						# print(f"delta_mu_init = {mu_equations([pu[0], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						if abs(p1[0]-p2[1])>1e-3 or abs(p1[1]-p2[0])>1e-3:
							print(f"There is a problem in the second half of secondary...", flush=True)
							print(f"seed was phi_p = {pu[1]}, guess was phi_s = {pu[0]}, phi = {pl[0], pl[1]}", flush=True)
							print(f"delta_mu_init = {mu_equations([pu[0], pl[0], pl[1]])}, delta_mu_fin = {mu_equations([root[0], root[1], root[2]])}")
						break						
		return sol_bin_up, sol_bin_down
	######################################################

	def partition_along_tangent(self, binodal_points, crit_point):
		tslope     = tangent.tangent2(self.vs, self.vc, self.vp, crit_point[0], crit_point[1], self.chi_pc, self.chi_ps, self.chi_sc, self.root_up_s, self.root_lo_s)
		tdir       = np.array([1, tslope])/np.sqrt(1+tslope**2)
		d_crit     = binodal_points - crit_point 
		crit_point = crit_point.reshape(-1,1)
		cross_prod = np.cross (d_crit, tdir)
		positive_side = cross_prod>0
		negative_side = cross_prod<0

		# pick the larger side 
		if np.sum(positive_side) > np.sum(negative_side):
			binodal_points = binodal_points[positive_side]
		elif np.sum(positive_side) < np.sum(negative_side):
			binodal_points = binodal_points[negative_side]
		else:
			print("Something's off... not changing anything.")

		return binodal_points
	######################################################

	def perform_secondary_searches(self, sol_bin_pos, sol_bin_neg, center, crit_point, p_gen):
		max_idx = np.argmax(np.linalg.norm(sol_bin_pos - crit_point, axis=1))
		positive_arm_tip = sol_bin_pos[max_idx]
		max_idx = np.argmax(np.linalg.norm(sol_bin_neg - crit_point, axis=1))
		negative_arm_tip = sol_bin_neg[max_idx]
		print(f"positive_arm_tip = {positive_arm_tip}", flush=True)
		print(f"negative_arm_tip = {negative_arm_tip}", flush=True)
		phi_positive, phi_negative, min_distances = B.perform_secondary_sweeps(positive_arm_tip, negative_arm_tip, center, crit_point, p_gen)
		print(f"Concluded secondary sweeps!", flush=True)

		# only keep fractions that are sufficiently close
		delta = 1 
		mask  = min_distances < delta
		phi_positive = phi_positive[mask]
		phi_negative = phi_negative[mask]

		# sort phi_positive and phi_negative by distance from crit_point
		d_from_crit = np.linalg.norm(phi_positive-crit_point, axis=1)
		sort_indices = np.argsort(d_from_crit)
		phi_positive = phi_positive[sort_indices]

		d_from_crit = np.linalg.norm(phi_negative-crit_point, axis=1)
		sort_indices = np.argsort(d_from_crit)
		phi_negative = phi_negative[sort_indices]

		# start performing binodal blasts
		print(f"Performing binodal blasts...", flush=True)
		_sol_bin_pos, _sol_bin_neg = B.perform_secondary_binodal_blast(phi_positive, phi_negative, center, crit_point)
		print("Completed secondary blast!", flush=True)
		sol_bin_pos = np.vstack((sol_bin_pos, _sol_bin_pos))
		sol_bin_neg = np.vstack((sol_bin_neg, _sol_bin_neg))	
		return sol_bin_pos, sol_bin_neg
	######################################################

	def build_bridges(self, sol_bin_pos, sol_bin_neg, center, crit_point, nadd):
		central_axis = (crit_point-center)/np.linalg.norm(crit_point-center)
		sol_bin_pos_loaded, m_load = ternary.add_rows_between_largest_gap(sol_bin_pos, nadd)
		sol_bin_neg_loaded         = ternary.add_rows_at_index(sol_bin_neg, m_load, nadd) 

		print(f"m_load = {m_load}")
		print(f"sol_bin_pos.shape = {sol_bin_pos.shape} and sol_bin_neg.shape = {sol_bin_neg.shape}",flush=True)
		print(f"sol_bin_pos_loaded.shape = {sol_bin_pos_loaded.shape} and sol_bin_neg.shape = {sol_bin_neg_loaded.shape}",flush=True)
		print(f"from {sol_bin_pos_loaded[m_load+1]} to {sol_bin_pos_loaded[m_load+nadd-1]}", flush=True)
		
		add_counter = 0
		for idx, pu in enumerate(sol_bin_pos_loaded[m_load+1:m_load+nadd-1]):
			def mu_equations(phi):
				eq1 = (self.mu_a(phi[0], pu[1])-self.mu_a(phi[1], phi[2])) # /(np.linalg.norm(np.array([phi[0], pu[1]])-np.array([phi[1],phi[2]])))
				eq2 = (self.mu_b(phi[0], pu[1])-self.mu_b(phi[1], phi[2])) # /(np.linalg.norm(np.array([phi[0], pu[1]])-np.array([phi[1],phi[2]])))
				eq3 = (self.mu_c(phi[0], pu[1])-self.mu_c(phi[1], phi[2])) # /(np.linalg.norm(np.array([phi[0], pu[1]])-np.array([phi[1],phi[2]])))
				return [eq1, eq2, eq3]
			
			root_store = []
			dist_store = []

			for jdx, pl in enumerate(sol_bin_neg_loaded[m_load+idx-nadd:m_load+1+idx+nadd]):
				root = fsolve(mu_equations, [pu[0], pl[0], pl[1]])
				if (np.abs(np.array(mu_equations(root))) > 1e-6).any():
					continue
				else:
					p1 = np.array([root[0], pu[1]])
					p2 = np.array([root[1], root[2]])
					if np.isnan(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isnan(ternary.stab_crit(p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue
					elif np.isinf(ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)) or np.isinf(ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)):
						continue
					elif args.ons and (ternary.stab_crit(p1[0], p1[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc)<0 or ternary.stab_crit (p2[0], p2[1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0):
						continue
					else:
						if np.sign(np.cross (central_axis, p1-center)) == np.sign(np.cross (central_axis, p2-center)):
							continue
						elif np.cross (central_axis, p1-center)>=0:
							root_store.append((p1,p2))
						else:
							root_store.append((p2,p1))
						dist_store.append(np.linalg.norm(p1-p2))		
			# choose the root that is furthest away
			try:
				best_root   = np.argmax(dist_store)
				root_combo  = root_store[best_root]
				sol_bin_pos = np.insert(sol_bin_pos, m_load+1+add_counter, root_combo[0], axis=0)
				sol_bin_neg = np.insert(sol_bin_neg, m_load+1+add_counter, root_combo[1], axis=0)
			except:
				print(f"Problem with a particular value. No solution found for pu = {pu}", flush=True)	

		return sol_bin_pos, sol_bin_neg	
	#!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!
	#						END OF CLASS
	#!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!

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

	# set up the simulation object
	B       = Binodal(inputs)
	B.obtain_crits()
	print(f"crits = {B.crits}", flush=True)
	# B.crits = np.array([[0.37037037, 0.37037037],[0.37037037,0.25925926],[0.25925926, 0.37037037]])
	B.stability_plots(ax)

	# perform the first sweep
	p_gen = dict()
	# DIVIDE BY 2 works alright!
	p_gen["rmin"]       = 0.01
	p_gen["rmax"]       = 0.1
	p_gen["r_dens"]     = 25
	p_gen["theta_dens"] = 25

	# define the center and crit point around which you want to find the binodal
	center     = np.mean(B.crits, axis=0)
	crit_point = B.crits[0]
	print(f"crit_point = {crit_point}", flush=True)

	phi_positive, phi_negative, min_distances = B.perform_first_sweep(center, crit_point, p_gen)
	# threshold min_distances
	delta = 1
	mask  = min_distances < delta
	phi_positive = phi_positive[mask]
	phi_negative = phi_negative[mask]

	# sort phi_positive and phi_negative by distance from crit_point
	d_from_crit = np.linalg.norm(phi_positive-crit_point, axis=1)
	sort_indices = np.argsort(d_from_crit)
	phi_positive = phi_positive[sort_indices]

	d_from_crit = np.linalg.norm(phi_negative-crit_point, axis=1)
	sort_indices = np.argsort(d_from_crit)
	phi_negative = phi_negative[sort_indices]

	# do the binodal blasts
	sol_bin_pos, sol_bin_neg = B.perform_first_binodal_blast(phi_positive, phi_negative, center, crit_point)

	# scale p_gen for secondary blasts
	p_gen["rmin"]       = 0.01/10
	p_gen["rmax"]       = 0.1/10

	# now start doing bigger iterations
	print(f"Start secondary blasts...", flush=True)
	for it in range(args.it):
		sol_bin_pos, sol_bin_neg = B.perform_secondary_searches(sol_bin_pos, sol_bin_neg, center, crit_point, p_gen)

	# sort phi_positive and phi_negative by distance from crit_point
	d_from_crit  = np.linalg.norm(sol_bin_pos-crit_point, axis=1)
	sort_indices = np.argsort(d_from_crit)
	sol_bin_pos  = sol_bin_pos[sort_indices]

	d_from_crit = np.linalg.norm(sol_bin_neg-crit_point, axis=1)
	sort_indices = np.argsort(d_from_crit)
	sol_bin_neg = sol_bin_neg[sort_indices]

	# now we can start filling up the gaps
	sol_bin_pos = np.vstack((crit_point, sol_bin_pos))
	sol_bin_neg = np.vstack((crit_point, sol_bin_neg))

	# build bridges from max_dist_idx to the next_idx
	print("Building bridges...", flush=True)
	nadd = 50
	for it in range(args.it):
		sol_bin_pos, sol_bin_neg = B.build_bridges(sol_bin_pos, sol_bin_neg, center, crit_point, nadd)

	print("Partitioning sol_bin_pos along tangent...", flush=True)
	sol_bin_pos = B.partition_along_tangent(sol_bin_pos, crit_point)
	print("Partitioning sol_bin_neg along tangent...", flush=True)
	sol_bin_neg = B.partition_along_tangent(sol_bin_neg, crit_point)

	ax.scatter(sol_bin_pos[:,0], 1-sol_bin_pos[:,0]-sol_bin_pos[:,1], sol_bin_pos[:,1], edgecolors='k', s=0.5, c='steelblue')
	ax.scatter(sol_bin_neg[:,0], 1-sol_bin_neg[:,0]-sol_bin_neg[:,1], sol_bin_neg[:,1], edgecolors='k', s=0.5, c='coral')

	# create the image
	if args.img != "None":
		if (".png" in args.img[-4:]):
			img_name = args.img
		elif ("." in args.img):
			img_name = args.img + ".png"
		else:
			img_name = args.img
		plt.savefig (img_name, dpi=1200, bbox_inches="tight")
	else:
		plt.savefig (f"bin_tern-vs_{B.vs}-vc_{B.vc}-vp_{B.vp}-chisc_{B.chi_sc}-chips_{B.chi_ps}-chipc_{B.chi_pc}.png", dpi=1200)

	print("Plotting gradient...", flush=True)
	fig_ = plt.figure(num=1, figsize=(3,3))
	ax_  = plt.axes()
	dy_dx_pos = np.gradient(sol_bin_pos[:,1])
	ax_.plot(sol_bin_pos[:,0], dy_dx_pos, c='darkred')
	dy_dx_neg = np.gradient(sol_bin_neg[:,1])
	ax_.plot(sol_bin_neg[:,0], dy_dx_neg, c='grey')
	fig_.savefig("gradient.png", dpi=1200, bbox_inches="tight")

	stop = time.time()
	print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)
