import numpy as np
import sympy as s

# rho_B is the number density of the bulk
# x_Cb is the cosolvent number density of the bulk
# rho_G is the number density of the gyration volume 
# rho_Gm is the number density of the monomer in the gyration volume
# rho_Gs is the number density of the solvent in the gyration volume
# rho_Gc is the number density of the cosolvent in the gyration volume
# rho_G = rho_Gm + rho_Gs + rho_Gc

class Bulk:
	def __init__(self, T, rho, x_c, sig_s, sig_c, e_ss, e_cc, e_sc):
		self.T     = T     # this is the bulk temperature
		self.rho   = rho   # this is the bulk number density 
		self.x_c   = x_c   # this is the cosolvent number fraction in the bulk
		self.sig_s = sig_s # this is the diameter of the hard sphere
		self.sig_c = sig_c # this is thed diameter of the hard sphere
		self.e_ss  = e_ss  # this is the SS interaction energy
		self.e_cc  = e_cc  # this is the CC interaction energy
		self.e_sc  = e_sc  # this is the SC interaction energy
		return 
	
	# define the exclusionary volumes 
	v_ss = lambda self: 4/3 * np.pi * self.sig_s ** 3
	v_cc = lambda self: 4/3 * np.pi * self.sig_c ** 3
	v_sc = lambda self: 4/3 * np.pi * ((self.sig_s + self.sig_c)/2) ** 3

	# define the attractive potentials 
	a_ss = lambda self: self.v_ss() * self.e_ss
	a_cc = lambda self: self.v_cc() * self.e_cc
	a_sc = lambda self: self.v_sc() * self.e_sc

	def xi_s(self):
		xi = np.pi * sig_s ** 3 / 6 * self.rho * (1 - self.x_c)
		return xi

	def xi_c(self):
		xi = np.pi * sig_c ** 3 / 6 * self.rho * (self.x_c)
		return xi

	def xi(self):
		return self.xi_s() + self.xi_c() 

	def delta_sc(self):
		delta = np.sqrt(self.xi_s() * self.xi_c())/self.xi() * (self.sig_c - self.sig_s) ** 2 /(self.sig_c * self.sig_s) * np.sqrt(self.x_c * (1-self.x_c))
		return delta

	def y1(self):
		return self.delta_sc() * (self.sig_s + self.sig_c)/np.sqrt(self.sig_s * self.sig_c)

	def y2(self):
		return 1/self.xi() * (self.xi_c()/self.sig_c + self.xi_s()/self.sig_s) * (self.delta_sc() * np.sqrt(self.sig_s * self.sig_c))

	def y3(self):
		return ((self.xi_c()/self.xi())**(2/3) * self.x_c**(1/3) + (self.xi_s()/self.xi())**(2/3) * (1-self.x_c)**(1/3))**3

	def A(self):
		term1 = -3/2*(1-self.y1() + self.y2() + self.y3())
		term2 = (3 * self.y2() + 2 * self.y3())/(1 - self.xi())
		term3 = 3/2 * (1 - self.y1() - self.y2() - self.y3()/3)/(1 - self.xi()) ** 2
		term4 = (self.y3() - 1) * np.log(1 - self.xi())
		return term1 + term2 + term3 + term4

	def pressure(self):
		term1 = (1 + self.xi() + self.xi()**2 - 3 * self.xi() * (self.y1() + self.y2() * self.xi() + self.xi()**2 * self.y3()/3))/(1-self.xi())**3
		term2 = self.rho/(2 * self.T) * (self.a_ss() * (1-self.x_c)**2 + self.a_cc() * self.x_c**2 + 2 * self.a_sc() * self.x_c * (1-self.x_c))
		return self.rho * self.T * (term1 - term2)

	def f_mix(self):
		term1 = self.rho * self.T * (self.x_c * (np.log(self.rho * self.x_c * 1/T**(3/2)) - 1) + (1 - self.x_c) * (np.log(self.rho * (1-self.x_c) * 1/T**(3/2)) - 1)) 
		term2 = self.rho * self.T * self.A()
		term3 = -1/2 * self.rho ** 2 *(self.a_ss() * (1 - self.x_c)**2 + self.a_cc() * x**2 + 2 * self.a_sc() * (1 - self.x_c) * self.x_c)
		return term1 + term2 + term3

	def df_mix(self):
		term1 = self.rho * self.T * (np.log(self.rho * self.x_c/self.T**(3/2))-np.log(self.rho * (1 - self.x_c)/self.T**(3/2)))
		term2 = -self.rho ** 2 * (self.a_sc() + self.a_ss() * (-1 + self.x_c) + self.a_cc() * self.x_c - 2 * self.a_sc() * self.x_c)
		return term1 + term2

	def mu_s(self):
		return 1/self.rho * (self.pressure() + self.f_mix() - self.x_c * self.df_mix())

	def mu_c(self):
		return 1/self.rho * (self.pressure() + self.f_mix() + (1 - self.x_c) * self.df_mix())


class Solvation:
	def __init__(self, T, P, alpha, Nm, Ns, Nc, sig_m, sig_s, sig_c, e_mm, e_ss, e_cc, e_ms, e_mc, e_sc):
		self.T     = T     # this is the bulk temperature
		self.P     = P     # this is the bulk pressure
		self.alpha = alpha # this is the Rg/Rg0 ratio
		self.Nm    = Nm    # number of monomers in the gyration volume
		self.Ns    = Ns    # number of solvent particles in the gyration volume
		self.Nc    = Nc    # number of cosolvent particles in the gyration volume
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
	V = lambda self: 2 * np.pi * self.Nm ** (3/2) * self.alpha ** 3 * /(9 * np.sqrt(6))

	# densities of the species
	rho_m = lambda self: self.Nm / self.V()
	rho_s = lambda self: self.Ns / self.V()
	rho_c = lambda self: self.Nc / self.V()
	rho   = lambda self: self.rho_m() + self.rho_s() + self.rho_c()

	def xi_m(self):
		xi = np.pi * sig_m ** 3 / 6 * self.rho_m()
		return xi

	def xi_s(self):
		xi = np.pi * sig_s ** 3 / 6 * self.rho_s()
		return xi

	def xi_c(self):
		xi = np.pi * sig_c ** 3 / 6 * self.rho_c()
		return xi

	def xi(self):
		return self.xi_m() + self.xi_s() + self.xi_c()

	def delta_ms(self):
		return np.sqrt(self.xi_m() * self.xi_s()) / self.xi() * (self.sig_s - self.sig_m)**2/(self.sig_s * self.sig_m) * np.sqrt(self.rho_s() * self.rho_m())/self.rho()

	def delta_mc(self):
		return np.sqrt(self.xi_m() * self.xi_c()) / self.xi() * (self.sig_c - self.sig_m)**2/(self.sig_c * self.sig_m) * np.sqrt(self.rho_c() * self.rho_m())/self.rho()

	def delta_sc(self):
		return np.sqrt(self.xi_s() * self.xi_c()) / self.xi() * (self.sig_s - self.sig_c)**2/(self.sig_s * self.sig_c) * np.sqrt(self.rho_s() * self.rho_c())/self.rho()

	def y1(self):
		term1 = self.delta_mc() * (self.sig_c + self.sig_m)/np.sqrt(self.sig_c * self.sig_m)
		term2 = self.delta_ms() * (self.sig_s + self.sig_m)/np.sqrt(self.sig_s * self.sig_m)
		term3 = self.delta_sc() * (self.sig_s + self.sig_c)/np.sqrt(self.sig_s * self.sig_c)
		return term1 + term2 + term3 

	def y2(self):
		term1 = 1/self.xi() * (self.xi_c()/self.sig_c + self.xi_s()/self.sig_s + self.xi_m()/self.sig_m) 
		term2 = (self.delta_mc() * np.sqrt(self.sig_m * self.sig_c) + self.delta_ms() * np.sqrt(self.sig_m * self.sig_s) + self.delta_sc() * np.sqrt(self.sig_s * self.sig_c))
		return term1 * term2

	def y3(self):
		term1 = (self.xi_c()/self.xi())**(2/3) * (self.rho_c()/self.rho())**(1/3)
		term2 = (self.xi_s()/self.xi())**(2/3) * (self.rho_s()/self.rho())**(1/3)
		term3 = (self.xi_m()/self.xi())**(2/3) * (self.rho_m()/self.rho())**(1/3)
		return (term1 + term2 + term3) ** 3

	def A(self):
		term1 = -3/2*(1-self.y1() + self.y2() + self.y3())
		term2 = (3 * self.y2() + 2 * self.y3())/(1 - self.xi())
		term3 = 3/2 * (1 - self.y1() - self.y2() - self.y3()/3)/(1 - self.xi()) ** 2
		term4 = (self.y3() - 1) * np.log(1 - self.xi())
		return term1 + term2 + term3 + term4

	def f_id(self):
		term1 = 9/4 * self.T * (self.alpha**2 + 1/(self.alpha**2))
		term2 = self.Nm * self.T * (np.log(self.rho_s() / self.T**(3/2)) - 1) + self.Ns * self.T * (np.log(self.rho_c()/self.T**(3/2)) - 1)
		return term1 + term2 

	def f_ex(self):
		term1 = self.rho() * self.T * self.A() 
		term2 = 1/2 * (self.a_mm() * self.rho_m()**2 + self.a_ss() * self.rho_s()**2 + self.a_cc() * self.rho_c()**2 + 2 * self.a_sc() * self.rho_s() * self.rho_c() + 2 * self.a_ms() * self.rho_m() * self.rho_s() + 2 * self.a_mc() * self.rho_m() * self.rho_c())
		return (term1 - term2) * self.V() 




if __name__=="__main__":

	# define the degree of polymerization
	Vg = lambda alpha, Nm: 2 * np.pi * Nm ** (3/2) * alpha ** 3 * b ** 3 /(9 * np.sqrt(6))

	# define the exclusionary volumes 
	v_ss = 4/3 * np.pi * sig_s ** 3
	v_cc = 4/3 * np.pi * sig_c ** 3
	v_mm = 4/3 * np.pi * sig_m ** 3
	v_sc = 4/3 * np.pi * ((sig_s + sig_c)/2) ** 3
	v_ms = 4/3 * np.pi * ((sig_m + sig_s)/2) ** 3
	v_mc = 4/3 * np.pi * ((sig_m + sig_c)/2) ** 3

	# define the attractive potentials 
	a_ss = v_ss * e_ss
	a_cc = v_cc * e_cc
	a_mm = v_mm * e_mm
	a_sc = v_sc * e_sc
	a_ms = v_ms * e_ms
	a_mc = v_mc * e_mc

	# define monomer density
	rho_i = lambda alpha, Ni, Nm: Ni / Vg(alpha, Nm)
	rho_G = lambda alpha, Nm, Ns, Nc: rho_i(alpha, Ns, Nm) + rho_i(alpha, Nc, Nm) + rho_i(alpha, Nm, Nm)

	# xi's are defined for the gyration volume
	# define a xi for each species, i
	xi_i   = lambda sig, rho: (np.pi * sig ** 3 * rho)/6
	xi_net = lambda alpha, Nm, Ns, Nc: xi_i(sig_s, rho_i(alpha, Ns, Nm)) + xi_i(sig_c, rho_i(alpha, Nc, Nm)) + xi_i(sig_m, rho_m(alpha, Nm, Nm))

	# delta's are also defined for the gyration volume 
	delta_ms = lambda alpha, Nm, Ns, Nc: np.sqrt(xi_i(sig_s, rho_i(alpha, Ns, Nm))*xi_i(sig_m, rho_i(alpha, Nm, Nm)))/xi_net(alpha, Nm, Ns, Nc) * (sig_s - sig_m)**2/(sig_s * sig_m) * np.sqrt(rho_i(alpha, Ns, Nm) * rho_i(alpha, Nm, Nm))/rho_G(alpha, Nm, Ns, Nc)
	delta_mc = lambda alpha, Nm, Ns, Nc: np.sqrt(xi_i(sig_c, rho_i(alpha, Nc, Nm))*xi_i(sig_m, rho_i(alpha, Nm, Nm)))/xi_net(alpha, Nm, Ns, Nc) * (sig_c - sig_m)**2/(sig_c * sig_m) * np.sqrt(rho_i(alpha, Nc, Nm) * rho_i(alpha, Nm, Nm))/rho_G(alpha, Nm, Ns, Nc)
	delta_sc = lambda alpha, Nm, Ns, Nc: np.sqrt(xi_i(sig_s, rho_i(alpha, Ns, Nm))*xi_i(sig_c, rho_i(alpha, Nc, Nm)))/xi_net(alpha, Nm, Ns, Nc) * (sig_s - sig_c)**2/(sig_s * sig_c) * np.sqrt(rho_i(alpha, Ns, Nm) * rho_i(alpha, Nc, Nm))/rho_G(alpha, Nm, Ns, Nc)
	
	y1 = lambda alpha, Nm, Ns, Nc: delta_mc(alpha, Nm, Ns, Nc) * (sig_c + sig_m)/np.sqrt(sig_c * sig_m) + delta_ms(alpha, Nm, Ns, Nc) * (sig_s + sig_m)/np.sqrt(sig_s*sig_s) + delta_sc(alpha, Nm, Ns, Nc) * (sig_s + sig_c)/np.sqrt(sig_s*sig_c)
	y2 = lambda alpha, Nm, Ns, Nc: 1/xi_net(alpha, Nm, Ns, Nc) * (xi_i(sig_c, rho_i(alpha, Nc, Nm))/sig_c + xi_i(sig_s, Ns, Nm)/sig_s + xi_m(alpha, Nm, Nm)/sig_m) * (delta_mc(alpha, Nm, Ns, Nc) * np.sqrt(sig_c * sig_m) + delta_ms(alpha, Nm, Ns, Nc) * np.sqrt(sig_s * sig_m) + delta_sc(alpha, Nm, Ns, Nc) * np.sqrt(sig_s, sig_c))
	y3 = lambda alpha, Nm, Ns, Nc: ((xi_i(sig_c, rho_i(sig_c, Nc, Nm))/xi_net(alpha, Nm, Ns, Nc))**2/3 * (rho_i(alpha, Nc, Nm)/rho_G(alpha, Nm, Ns, Nc))**1/3 + (xi_i(sig_s, rho_i(sig_s, Ns, Nm))/xi_net(alpha, Nm, Ns, Nc))**2/3 * (rho_i(alpha, Ns, Nm)/rho_G(alpha, Nm, Ns, Nc))**1/3 + (xi_i(sig_m, rho_i(sig_m, Nm, Nm))/xi_net(alpha, Nm, Ns, Nc))**2/3 * (rho_i(alpha, Nm, Nm)/rho_G(alpha, Nm, Ns, Nc))**1/3) ** 3

	A = lambda alpha, Nm, Ns, Nc: -3/2 * (1 - y1(alpha, Nm, Ns, Nc) + y2(alpha, Nm, Ns, Nc) + y3(alpha, Nm, Ns, Nc)) + (3 * y2 (alpha, Nm, Ns, Nc) + 2 * y3(alpha, Nm, Ns, Nc))/(1 - xi_net(alpha, Nm, Ns, Nc)) \
	+ 3 * (1-y1(alpha, Nm, Ns, Nc) - y2(alpha, Nm, Ns, Nc) - y3(alpha, Nm, Ns, Nc)/3)/(2 * (1-xi_net(alpha, Nm, Ns, Nc))) ** 2 + (y3(alpha, Nm, Ns, Nc) -1) * np.log(1 - xi_net(alpha, Nm, Ns, Nc))

	pressure = lambda rho_B, T, x: 
