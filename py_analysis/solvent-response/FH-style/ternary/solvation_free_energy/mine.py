import numpy as np
import cma

# rho_B is the number density of the bulk
# x_Cb is the cosolvent number density of the bulk
# rho_G is the number density of the gyration volume 
# rho_Gm is the number density of the monomer in the gyration volume
# rho_Gs is the number density of the solvent in the gyration volume
# rho_Gc is the number density of the cosolvent in the gyration volume
# rho_G = rho_Gm + rho_Gs + rho_Gc

z = 26

class Bulk:
	def __init__(self, T, x_c, e_ss, e_cc, e_sc):
		self.T     = T     # this is the bulk temperature
		self.x_c   = x_c   # this is the cosolvent number fraction in the bulk
		self.e_ss  = e_ss  # this is the SS interaction energy
		self.e_cc  = e_cc  # this is the CC interaction energy
		self.e_sc  = e_sc  # this is the SC interaction energy
		self.chi   = (z-2)/T * (e_sc - 1/2 * (e_ss + e_cc))
		return 

	def mu_s(self):
		mu = np.log(1-self.x_c) + self.chi * (self.x_c)**2
		return mu

	def mu_c(self):
		mu = np.log(self.x_c) + self.chi * (1-self.x_c)**2
		return mu

class Gyration:
	def __init__(self, T, Nm, sig_m, sig_s, sig_c, e_mm, e_ss, e_cc, e_ms, e_mc, e_sc):
		self.T     = T     # this is the bulk temperature
		self.Nm    = Nm    # number of monomers in the gyration volume
		self.sig_m = sig_m # this is the diameter of the hard sphere
		self.sig_s = sig_s # this is the diameter of the hard sphere
		self.sig_c = sig_c # this is the diameter of the hard sphere
		self.e_mm  = e_mm  # this is the MM interaction energy
		self.e_ss  = e_ss  # this is the SS interaction energy
		self.e_cc  = e_cc  # this is the CC interaction energy
		self.e_ms  = e_ms  # this is the MS interaction energy
		self.e_mc  = e_mc  # this is the MC interaction energy
		self.e_sc  = e_sc  # this is the SC interaction energy
		self.chi_sc = (z-2)/T * (e_sc - 1/2 * (e_ss + e_cc))
		self.chi_ms = (z-2)/T * (e_ms - 1/2 * (e_mm + e_ss))
		self.chi_mc = (z-2)/T * (e_mc - 1/2 * (e_mm + e_cc))
		return 

	# define the volume
	V = lambda self, alpha: 2 * np.pi * (self.Nm ** (3/2)) * alpha ** 3/(9 * np.sqrt(6))

	def f_id(self, alpha):
		term1 = 9/4 * self.T * (alpha**2 + 1/(alpha**2))
		return term1

	def f_ex(self, alpha, Ns, Nc):
		Vs = Ns * np.pi/6 * (self.sig_s)**3
		Vc = Nc * np.pi/6 * (self.sig_c)**3
		Vm = Nm * np.pi/6 * (self.sig_m)**3
		entropy = Ns * np.log(Vs/self.V(alpha)) + Nc * np.log(Vc/self.V(alpha)) + Nm * np.log(Vm/self.V(alpha))
		energy  = self.chi_sc * Ns * Vs/self.V(alpha) + self.chi_mc * self.Nm * Vc/self.V(alpha) + self.chi_ms * self.Nm * Vs/self.V(alpha)
		return self.T * (entropy + energy)

if __name__=="__main__":

	# initialize the bulk variables
	T     = 0.1
	x_c   = 0.5
	e_ss  = 0
	e_cc  = 0
	e_sc  = -2
	bulk = Bulk(T, x_c, e_ss, e_cc, e_sc)
	
	# initialize the variables in the gyration volume
	Nm    = 50
	sig_m = 1
	sig_s = 1
	sig_c = 1
	e_mm  = 0
	e_ms  = 0
	e_mc  = 0
	gyr = Gyration(T, Nm, sig_m, sig_s, sig_c, e_mm, e_ss, e_cc, e_ms, e_mc, e_sc)

	print(f"mu_s = {bulk.mu_s()}")
	print(f"mu_c = {bulk.mu_c()}")

	def Gsolv(params):
		if params[0] < 1e-6 or params[0] > 4 or params[1] < 1e-6 or params[2] < 1e-6 or params[3] > 1e-6:
			obj = 1e+8
		else:
			G = gyr.f_id(params[0]) + gyr.f_ex(params[0], params[1], params[2]) - T * params[1] * bulk.mu_s() - T * params[2] * bulk.mu_c() + params[3] * gyr.V(params[0])
			if np.isnan(G) or np.isinf(G):
				obj = 1e+8
			else:
				obj = G
		return obj
	
	x0 = [1, 1, 1, 1]
	result = cma.fmin(Gsolv, x0, sigma0=0.01)

	# print the result
	print(f"Minimum value: {result[1]}",   flush=True)
	print(f"Optimal solution: {result[0]}", flush=True)
	