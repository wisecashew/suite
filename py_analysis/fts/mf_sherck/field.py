#!/home/satyend/.conda/envs/phase/bin/python
import numpy as np
kB = 1

class Field:
	def __init__(self, inputs):
		self.v_pp     = inputs["v_pp"]
		self.v_pa     = inputs["v_pa"]
		self.v_pb     = inputs["v_pb"]
		self.v_aa     = inputs["v_aa"]
		self.v_bb     = inputs["v_bb"]
		self.v_ab     = inputs["v_ab"]
		self.a_p      = inputs["a_p" ]
		self.a_a      = inputs["a_a" ]
		self.a_b      = inputs["a_b" ]
		self.deg_poly = inputs["deg_poly"]
		self.rho0     = inputs["rho0"]
		self.T        = np.logspace(-6, 2, 500)
		self.beta     = 1 / (kB * self.T)
		self.get_self_interactions()
		self.get_cross_interactions()
		return
	
	def print(self):
		print(f"v_pp: {self.v_pp}")
		print(f"v_pa: {self.v_pa}")
		print(f"v_pb: {self.v_pb}")
		print(f"v_aa: {self.v_aa}")
		print(f"v_bb: {self.v_bb}")
		print(f"v_ab: {self.v_ab}")
		print(f"a_p : {self.a_p }")
		print(f"a_a : {self.a_a }")
		print(f"a_b : {self.a_b }")
		print(f"v_ps: {self.v_ps}")
		print(f"v_ss: {self.v_ss}")
		print(f"deg_poly: {self.deg_poly}")
		return

	def get_self_interactions(self):
		self.V_pp0 = self.v_pp / (4 * np.pi * self.a_p**2)**(3/2)
		self.V_aa0 = self.v_aa / (4 * np.pi * self.a_a**2)**(3/2)
		self.V_bb0 = self.v_bb / (4 * np.pi * self.a_b**2)**(3/2)
		return
	
	def get_cross_interactions(self):
		self.v_ps = self.v_pa + self.v_pb
		self.v_ss = 1/2 * (self.v_aa + self.v_bb + 2 * self.v_ab)
		return

	def chem_pot_P(self, rho_p, rho_s):
		# this computes the chemical potential of P in a phase
		# it ignores the constants present in the expression
		entr_p            = np.log(rho_p/self.deg_poly)
		self_intr         = - 1 / 2 * self.beta * self.deg_poly * self.V_pp0
		other_interaction = self.beta * self.deg_poly * (self.v_pa * rho_s + self.v_pb * rho_s + self.v_pp * rho_p)
		mu_P              = entr_p + self_intr + other_interaction
		return mu_P

	def chem_pot_S(self, rho_p, rho_s):
		# this computes the chemical potential of S in a phase
		# it ignores the constants present in the expression
		entr_s            = np.log(rho_s)
		self_intr         = - 1 / 2 * self.beta * (self.V_aa0 + self.V_bb0)
		other_interaction = self.beta * (self.v_ps * rho_p + 2 * self.v_ss * rho_s)
		mu_S              = entr_s + self_intr + other_interaction
		return mu_S

	def stability(self, rho_p, rho_s):
		d2f_dp2  = 1/(rho_p * self.deg_poly) + 1 * self.beta * self.v_pp
		d2f_ds2  = 1/rho_s                   + 1 * self.beta * self.v_ss
		d2f_dpds = self.beta * self.v_ps
		stab     = d2f_dp2 * d2f_ds2 - d2f_dpds**2
		return stab

	def stab_roots(self, rho_p, rho_s):
		# print(f"Leading term is ", end='', flush=True)
		# if 2 * self.v_pp * self.v_ss - self.v_ps**2 > 0:
		#	print(f"positive.", flush=True)
		# else:
		#	print(f"negative.", flush=True)
		
		# moving on...
		roots = []
		discriminant = (self.v_ss/(rho_p * self.deg_poly) + self.v_pp/rho_s)**2 - \
			4 *(self.v_pp * self.v_ss - self.v_ps**2) * 1 / (rho_p * rho_s * self.deg_poly)

		# print(f"Discriminant is {discriminant}", flush=True)
		if discriminant < 0 or np.isnan(discriminant or np.isinf(discriminant)):
			print(f"No temperatures found.")
			return roots
		else:
			beta1 = (-(self.v_ss/(rho_p * self.deg_poly) + self.v_pp/rho_s) - np.sqrt(discriminant)) / \
				(2 * (self.v_pp * self.v_ss - self.v_ps**2))
			if beta1 > 0:
				roots.append(beta1)
			#################
			beta2 = (-(self.v_ss/(rho_p * self.deg_poly) + self.v_pp/rho_s) + np.sqrt(discriminant)) / \
				(2 * (self.v_pp * self.v_ss - self.v_ps**2))
			if beta2 > 0:
				roots.append(beta2)
			return roots

