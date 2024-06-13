import numpy as np
from scipy.optimize import minimize

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
		# self.Ns    = Ns    # number of solvent particles in the gyration volume
		# self.Nc    = Nc    # number of cosolvent particles in the gyration volume
		self.sig_m = sig_m # this is the diameter of the hard sphere
		self.sig_s = sig_s # this is the diameter of the hard sphere
		self.sig_c = sig_c # this is the diameter of the hard sphere
		self.e_mm  = e_mm  # this is the MM interaction energy
		self.e_ss  = e_ss  # this is the SS interaction energy
		self.e_cc  = e_cc  # this is the CC interaction energy
		self.e_ms  = e_ms  # this is the MS interaction energy
		self.e_mc  = e_mc  # this is the MC interaction energy
		self.e_sc  = e_sc  # this is the SC interaction energy
		return 

	# define the exclusionary volumes
	v_mm = lambda self: 4/3 * np.pi * self.sig_m ** 3
	v_ss = lambda self: 4/3 * np.pi * self.sig_s ** 3
	v_cc = lambda self: 4/3 * np.pi * self.sig_c ** 3
	v_ms = lambda self: 4/3 * np.pi * ((self.sig_m + self.sig_s)/2) ** 3
	v_mc = lambda self: 4/3 * np.pi * ((self.sig_m + self.sig_c)/2) ** 3
	v_sc = lambda self: 4/3 * np.pi * ((self.sig_s + self.sig_c)/2) ** 3

	# define the attractive potentials 
	a_mm = lambda self: self.v_ss() * self.e_mm
	a_ss = lambda self: self.v_ss() * self.e_ss
	a_cc = lambda self: self.v_cc() * self.e_cc
	a_ms = lambda self: self.v_ms() * self.e_ms
	a_mc = lambda self: self.v_mc() * self.e_mc
	a_sc = lambda self: self.v_sc() * self.e_sc

	# define the volume
	V = lambda self, alpha: 2 * np.pi * (self.Nm ** (3/2)) * alpha ** 3/(9 * np.sqrt(6))

	# densities of the species
	rho_m = lambda self, alpha        : self.Nm / self.V(alpha)
	rho_s = lambda self, alpha, Ns    : Ns / self.V(alpha)
	rho_c = lambda self, alpha, Nc    : Nc / self.V(alpha)
	rho   = lambda self, alpha, Ns, Nc: self.rho_m(alpha) + self.rho_s(alpha, Ns) + self.rho_c(alpha, Nc)

	def xi_m(self, alpha):
		xi = np.pi * self.sig_m ** 3 / 6 * self.rho_m(alpha)
		return xi

	def xi_s(self, alpha, Ns):
		xi = np.pi * self.sig_s ** 3 / 6 * self.rho_s(alpha, Ns)
		return xi

	def xi_c(self, alpha, Nc):
		xi = np.pi * self.sig_c ** 3 / 6 * self.rho_c(alpha, Nc)
		return xi

	def xi(self, alpha, Ns, Nc):
		return self.xi_m(alpha) + self.xi_s(alpha, Ns) + self.xi_c(alpha, Nc)

	def delta_ms(self, alpha, Ns, Nc):
		return np.sqrt(self.xi_m(alpha) * self.xi_s(alpha, Ns)) / self.xi(alpha, Ns, Nc) * (self.sig_s - self.sig_m)**2/(self.sig_s * self.sig_m) * np.sqrt(self.rho_s(alpha, Ns) * self.rho_m(alpha))/self.rho(alpha, Ns, Nc)

	def delta_mc(self, alpha, Ns, Nc):
		return np.sqrt(self.xi_m(alpha) * self.xi_c(alpha, Nc)) / self.xi(alpha, Ns, Nc) * (self.sig_c - self.sig_m)**2/(self.sig_c * self.sig_m) * np.sqrt(self.rho_c(alpha, Nc) * self.rho_m(alpha))/self.rho(alpha, Ns, Nc)

	def delta_sc(self, alpha, Ns, Nc):
		return np.sqrt(self.xi_s(alpha, Ns) * self.xi_c(alpha, Nc)) / self.xi(alpha, Ns, Nc) * (self.sig_s - self.sig_c)**2/(self.sig_s * self.sig_c) * np.sqrt(self.rho_s(alpha, Ns) * self.rho_c(alpha, Nc))/self.rho(alpha, Ns, Nc)

	def y1(self, alpha, Ns, Nc):
		term1 = self.delta_mc(alpha, Ns, Nc) * (self.sig_c + self.sig_m)/np.sqrt(self.sig_c * self.sig_m)
		term2 = self.delta_ms(alpha, Ns, Nc) * (self.sig_s + self.sig_m)/np.sqrt(self.sig_s * self.sig_m)
		term3 = self.delta_sc(alpha, Ns, Nc) * (self.sig_s + self.sig_c)/np.sqrt(self.sig_s * self.sig_c)
		return term1 + term2 + term3 

	def y2(self, alpha, Ns, Nc):
		term1 = 1/self.xi(alpha, Ns, Nc) * (self.xi_c(alpha, Nc)/self.sig_c + self.xi_s(alpha, Ns)/self.sig_s + self.xi_m(alpha)/self.sig_m) 
		term2 = (self.delta_mc(alpha, Ns, Nc) * np.sqrt(self.sig_m * self.sig_c) + self.delta_ms(alpha, Ns, Nc) * np.sqrt(self.sig_m * self.sig_s) + self.delta_sc(alpha, Ns, Nc) * np.sqrt(self.sig_s * self.sig_c))
		return term1 * term2

	def y3(self, alpha, Ns, Nc):
		term1 = (self.xi_c(alpha, Nc)/self.xi(alpha, Ns, Nc))**(2/3) * (self.rho_c(alpha, Nc)/self.rho(alpha, Ns, Nc))**(1/3)
		term2 = (self.xi_s(alpha, Ns)/self.xi(alpha, Ns, Nc))**(2/3) * (self.rho_s(alpha, Ns)/self.rho(alpha, Ns, Nc))**(1/3)
		term3 = (self.xi_m(alpha)/self.xi(alpha, Ns, Nc))**(2/3) * (self.rho_m(alpha)/self.rho(alpha, Ns, Nc))**(1/3)
		return (term1 + term2 + term3) ** 3

	def A(self, alpha, Ns, Nc):
		term1 = -3/2*(1-self.y1(alpha, Ns, Nc) + self.y2(alpha, Ns, Nc) + self.y3(alpha, Ns, Nc))
		term2 = (3 * self.y2(alpha, Ns, Nc) + 2 * self.y3(alpha, Ns, Nc))/(1 - self.xi(alpha, Ns, Nc))
		term3 = 3/2 * (1 - self.y1(alpha, Ns, Nc) - self.y2(alpha, Ns, Nc) - self.y3(alpha, Ns, Nc)/3)/(1 - self.xi(alpha, Ns, Nc)) ** 2
		term4 = (self.y3(alpha, Ns, Nc) - 1) * np.log(1 - self.xi(alpha, Ns, Nc))
		return term1 + term2 + term3 + term4

	def f_id(self, alpha, Ns, Nc):
		term1 = 9/4 * self.T * (alpha**2 + 1/(alpha**2))
		term2 = Ns * self.T * (np.log(self.rho_s(alpha, Ns) / self.T**(3/2)) - 1) + Nc * self.T * (np.log(self.rho_c(alpha, Nc)/self.T**(3/2)) - 1)
		return term1 + term2 

	def f_ex(self, alpha, Ns, Nc):
		term1 = self.rho(alpha, Ns, Nc) * self.T * self.A(alpha, Ns, Nc) 
		term2 = 1/2 * (self.a_mm() * self.rho_m(alpha)**2 + \
				self.a_ss() * self.rho_s(alpha, Ns)**2 + \
				self.a_cc() * self.rho_c(alpha, Nc)**2 + \
				2 * self.a_sc() * self.rho_s(alpha, Ns) * self.rho_c(alpha, Nc) + \
				2 * self.a_ms() * self.rho_m(alpha) * self.rho_s(alpha, Ns) + \
				2 * self.a_mc() * self.rho_m(alpha) * self.rho_c(alpha, Nc))
		return (term1 - term2) * self.V(alpha) 

if __name__=="__main__":

	# initialize the bulk variables
	T     = 0.1
	x_c   = 0.5
	e_ss  = 0
	e_cc  = 0
	e_sc  = 2
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
	
	def con_size(params):
		return np.floor((params[0]/5))
	
	def con_positive(params):
		return ((params[0]>0) and (params[1]>0) and (params[2]>0))-1
	
	
	cons = [{"type":"eq", "fun": con_size}, {"type":"eq", "fun": con_positive}]

	def Gsolv(params):
		G = gyr.f_id(params[0], params[1], params[2]) + gyr.f_ex(params[0], params[1], params[2]) - params[1] * bulk.mu_s() - params[2] * bulk.mu_c() 
		return G
	
	x0 = [2, 1, 1]
	result = minimize(Gsolv, x0, constraints=cons)

	# print the result
	print(f"Minimum value: {result.fun}",   flush=True)
	print(f"Optimal solution: {result.x}", flush=True)
	