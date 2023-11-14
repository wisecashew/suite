import numpy as np
import tangent
from scipy.optimize import fsolve

def add_rows_between_largest_gap (array, M):

	diff        = np.diff (array, axis=0)
	distances   = np.linalg.norm (diff, axis=1)
	max_dists   = np.argmax (distances)
	insert_rows = np.linspace (array[max_dists], array[max_dists+1], M)
	array       = np.insert (array, max_dists+1, insert_rows[1:-1], axis=0)

	return array, max_dists

def add_rows_at_index (array, idx, M):

	insert_rows = np.linspace (array[idx], array[idx+1], M)
	array       = np.insert (array, idx+1, insert_rows[1:-1],axis=0)

	return array


def remove_close_rows(array, threshold):

	filtered_array = np.empty ((0,2))
	for i, elem in enumerate(array):
		if i == 0:
			filtered_array = np.vstack((filtered_array, elem))
			continue
		else:
			sieve = (np.linalg.norm(filtered_array - elem, axis=1) < 1e-3).any()
			if sieve:
				continue
			else:
				filtered_array = np.vstack((filtered_array, elem))

	return filtered_array


def crit_condition (vs, vc, vp, phi_p, phi_s, chi_sc, chi_ps, chi_pc):

	phi_c = 1-phi_p-phi_s
	t1    = 1/(vc*phi_c) + 1/(phi_p*vp) - 2*chi_pc
	t2    = (1/(vc*(phi_c)**2) - 1/(vs*phi_s**2))*(1/(vc*(phi_c)) + 1/(phi_p*vp) - 2*chi_pc) + (1/(vc*phi_c) + 1/(vs*phi_s) - 2*chi_sc)/(vc*(phi_c)**2) - 2*(1/(vc*phi_c) - chi_pc - chi_sc + chi_ps)/(vc*(phi_c)**2)

	u1    = (1/(vc*phi_c) + 1/(phi_p*vp) - 2*chi_pc)/(vc*(phi_c)**2) + (1/(vc*phi_c**2) - 1/(phi_p**2 * vp))*(1/(vc*phi_c) + 1/(vs*phi_s) - 2*chi_sc) - 2*(1/(vc*phi_c) + chi_ps - chi_sc - chi_pc)/(vc*phi_c**2)
	u2    = 1/(vc*phi_c) - chi_pc - chi_sc + chi_ps

	return t1*t2 - u1*u2

def find_crit_point (vs, vc, vp, chi_sc, chi_ps, chi_pc, root_up_p, root_up_s, root_lo_p, root_lo_s):

	def send_to_fsolve_r1 (phi_s):
		phi_p_upper = root_up_s (phi_s)
		return crit_condition (vs, vc, vp, phi_p_upper, phi_s, chi_sc, chi_ps, chi_pc)

	def send_to_fsolve_r2 (phi_s):
		phi_p_lower = root_lo_s (phi_s)
		return crit_condition (vs, vc, vp, phi_p_lower, phi_s, chi_sc, chi_ps, chi_pc)

	def send_to_fsolve_r3 (phi_p):
		phi_s_upper = root_up_p (phi_p)
		return crit_condition (vs, vc, vp, phi_p, phi_s_upper, chi_sc, chi_ps, chi_pc)

	def send_to_fsolve_r4 (phi_p):
		phi_s_lower = root_lo_p (phi_p)
		return crit_condition (vs, vc, vp, phi_p, phi_s_lower, chi_sc, chi_ps, chi_pc)

	guesses = np.linspace (0, 1, 10000)
	roots_up   = np.empty ((0,2))
	roots_down = np.empty ((0,2))


	print ("\tIn r1...", flush=True)
	for g in guesses:
		root = fsolve (send_to_fsolve_r1, g)

		if abs(send_to_fsolve_r1(root)) < 1e-12:

			if root >= 1 or root <= 0 or np.isnan(root):
				pass
			else:
				r_up  = root_up_s (root)[0]
				r_tup = np.array([root[0], r_up])
				if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_up):
					pass

				elif r_tup in roots_up:
					pass

				else:
					if len(roots_up) == 0:
						roots_up = np.vstack ((roots_up,r_tup))
					else:
						similarity = (np.linalg.norm(roots_up - r_tup, axis=1) < 1e-3).any()
						if similarity:
							pass
						else:
							roots_up = np.vstack ((roots_up,r_tup))   
		else:
			pass

	print ("\tIn r3...", flush=True)
	for g in guesses:
		root = fsolve (send_to_fsolve_r3, g)

		if abs(send_to_fsolve_r3(root)) < 1e-12:

			if root >= 1 or root <= 0 or np.isnan(root):
				pass
			else:
				r_up  = root_up_p (root)[0]
				r_tup = np.array([r_up, root[0]])
				if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_up):
					pass

				elif r_tup in roots_up:
					pass

				else:
					if len(roots_up) == 0:
						roots_up = np.vstack ((roots_up,r_tup))
					else:
						similarity = (np.linalg.norm(roots_up - r_tup, axis=1) < 1e-3).any()
						if similarity:
							pass
						else:
							roots_up = np.vstack ((roots_up,r_tup))   
		else:
			pass


	print ("\tIn r2...", flush=True)
	for g in guesses:
		root = fsolve (send_to_fsolve_r2, g)

		if abs(send_to_fsolve_r2(root)) < 1e-12:

			if root >= 1 or root <= 0 or np.isnan(root):
				pass
			else:
				r_lo = root_lo_s(root)[0]
				r_tup = np.array([root[0], r_lo])
				if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_lo):
					pass

				elif r_tup in roots_down:
					pass

				else: 
					if len(roots_down) == 0:
						roots_down = np.vstack ((roots_down,r_tup))
					else:
						similarity = (np.linalg.norm(roots_down - r_tup, axis=1) < 1e-3).any ()
						if similarity:
							pass
						else:
							roots_down = np.vstack ((roots_down,r_tup))

		else:
			pass


	print ("\tIn r4...", flush=True)
	for g in guesses:
		root = fsolve (send_to_fsolve_r4, g)

		if abs(send_to_fsolve_r4(root)) < 1e-12:

			if root >= 1 or root <= 0 or np.isnan(root):
				pass
			else:
				r_lo = root_lo_p(root)[0]
				r_tup = np.array([r_lo, root[0]])
				if (r_tup >= 1).any() or (r_tup <= 0).any() or np.sum(r_tup) >= 1 or np.isnan(r_lo):
					pass

				elif r_tup in roots_down:
					pass

				else:
					if len(roots_down) == 0:
						roots_down = np.vstack ((roots_down,r_tup))
					else:
						similarity = (np.linalg.norm(roots_down - r_tup, axis=1) < 1e-3).any ()
						if similarity:
							pass
						else:
							roots_down = np.vstack ((roots_down,r_tup))

		else:
			pass

	return roots_up, roots_down

#############################################################################################################

def stab_crit (p_s, p_p, vs, vc, vp, c_ps, c_pc, c_sc):
	return (1/(vp*p_p) + 1/(vc*(1-p_s - p_p)) - 2 * c_pc) * (1/(vs*p_s) + \
	1/(vc*(1-p_s - p_p)) - 2 * c_sc) - (1/(vc*(1-p_s-p_p)) + c_ps - c_pc - c_sc) ** 2

#############################################################################################################

def embelish(ax, tern_b):

	if tern_b:

		ax.set_tlabel ("$\\phi _{S}$") # ('Vol. frac. A')
		ax.set_llabel ("$\\phi _{C}$") # ('Vol. frac. C')
		ax.set_rlabel ("$\\phi _{P}$") # ('Vol. frac. B')
		ax.set_tlim(0, 1)
		ax.set_llim(0, 1)
		ax.set_rlim(0, 1)
		positions = ['tick1', 'tick2']
		for position in positions:
			ax.taxis.set_ticks_position(position)
			ax.laxis.set_ticks_position(position)
			ax.raxis.set_ticks_position(position)

	else:
		ax.set_xlabel ("$\\phi _{S}$")
		ax.set_ylabel ("$\\phi _{P}$")
		ax.set_xlim   (0,1)
		ax.set_ylim   (0,1)

	return

#############################################################################################################

def add_tang_norm(ax, tang_b, tern_b, crit_b, crits, vs, vc, vp, chi_pc, chi_ps, chi_sc, root_up_s, root_lo_s):

	if tang_b and crit_b:
		cols = ["steelblue", "limegreen", "pink", "maroon"]
		for idx,crit in enumerate(crits):    
			idx = idx%len(crits)
			slope               = tangent.tangent2 (vs, vc, vp, crit[0], crit[1], chi_pc, chi_ps, chi_sc, root_up_s, root_lo_s)
			perp_slope          = -1/slope
			tangent_vector      = np.array([1, slope]) / np.sqrt(1+slope**2)
			normal_vector       = np.array([1, perp_slope]) / np.sqrt(1+perp_slope**2)

			points_along_tangent = lambda L: np.array([crit[0] + L * tangent_vector[0], crit[1] + L * tangent_vector[1]])
			points_along_normal  = lambda L: np.array([crit[0] + L * normal_vector [0], crit[1] + L * normal_vector[1] ])

			l = np.linspace (-10, 10, 100)
			if tern_b:
				ax.plot (points_along_tangent(l)[0], 1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_tangent(l)[1], c=cols[idx], lw=1)
				ax.plot (points_along_normal(l)[0],  1-points_along_tangent(l)[0]-points_along_tangent(l)[1], points_along_normal(l)[1],  c=cols[idx], lw=1)

			else:
				ax.plot (points_along_tangent(l)[0],  points_along_tangent(l)[1],  c=cols[idx], lw=1)
				ax.plot (points_along_normal (l)[0],  points_along_normal(l) [1],  c=cols[idx], lw=1)

	return

#############################################################################################################

def plot(ax, tern_b, edges_b, crits_b, crits, chi_ps, chi_pc, chi_sc, p_s, p_p, cols, root_up_s, root_lo_s):

	if tern_b:
		ax.scatter(p_s, 1-p_p-p_s, p_p, s=1, color=cols)

		if crits_b:
			ax.scatter(crits[:,0], 1-crits[:,0]-crits[:,1], crits[:,1], color='darkred', edgecolors='darkred',s=4, zorder=15)
			ax.scatter(np.mean(crits, axis=0)[0], 1-np.mean(crits, axis=0)[0]-np.mean(crits, axis=0)[1], np.mean(crits, axis=0)[1], color="limegreen", edgecolors="limegreen", s=8, zorder=15)

		else:
			pass

		if edges_b:
			meshsize            = 1000
			phi_s               = np.linspace (0.001, 1-0.001, meshsize*10)
			chi_ps              = chi_ps
			chi_pc              = chi_pc
			chi_sc              = chi_sc

			r1 = root_up_s (p_s)
			r2 = root_lo_s (p_s)

			to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
			r1        = r1[to_keep_1]

			to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
			r2        = r2[to_keep_2]

			# Plot the points
			ax.scatter (p_s[to_keep_1], 1-p_s[to_keep_1]-r1, r1, color='springgreen',   s=1)
			ax.scatter (p_s[to_keep_2], 1-p_s[to_keep_2]-r2, r2, color='darkturquoise', s=1)

	else:
		ax.scatter (p_s, p_p, s=1, color=cols)
		if crits_b:
			ax.scatter (crits[:,0], crits[:,1], color='darkred', edgecolors='darkred', s=4, zorder=15)
			ax.scatter (np.mean (crits, axis=0)[0], np.mean(crits, axis=0)[1], color='limegreen', edgecolors='limegreen', s=8, zorder=15)
		else:
			pass

		if edges_b:
			meshsize            = 1000
			phi_s               = np.logspace (-6, np.log10(1-1e-6), meshsize*100)
			chi_ps              = chi_ps
			chi_pc              = chi_pc
			chi_sc              = chi_sc

			r1 = root_up_s (p_s, chi_ps, chi_pc, chi_sc)
			r2 = root_lo_s (p_s, chi_ps, chi_pc, chi_sc)

			to_keep_1 = (~np.isnan(r1)) * (r1 <= 1) * (r1 >= 0)
			r1 = r1[to_keep_1]

			to_keep_2 = (~np.isnan(r2)) * (r2 <= 1) * (r2 >= 0)
			r2 = r2[to_keep_2]

			# Plot the points
			ax.scatter(p_s[to_keep_1], r1, color='springgreen',   s=1)
			ax.scatter(p_s[to_keep_2], r2, color='darkturquoise', s=1)

	return

#############################################################################################################


