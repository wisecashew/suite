import numpy

class Mu_PS:
	def __init__(self, inputs):
		self.chi_sc = inputs["chi_sc"]
		self.chi_ps = inputs["chi_ps"]
		self.chi_pc = inputs["chi_pc"]
		self.vs     = inputs["vs"]
		self.vc     = inputs["vc"]
		self.vp     = inputs["vp"]
		return

	# calculate chemical potentials
	mu_s = lambda self, phi_s, phi_p: np.float64(np.log(phi_s)         + 1 - phi_s - self.vs/self.vp * phi_p - self.vs/self.vc * (1-phi_s-phi_p) + self.vs * (phi_p**2 * self.chi_ps + (1-phi_s-phi_p)**2 * self.chi_sc + phi_p * (1-phi_s-phi_p) * (self.chi_ps + self.chi_sc - self.chi_pc) ))
	mu_p = lambda self, phi_s, phi_p: np.float64(np.log(phi_p)         + 1 - phi_p - self.vp/self.vs * phi_s - self.vp/self.vc * (1-phi_s-phi_p) + self.vp * (phi_s**2 * self.chi_ps + (1-phi_s-phi_p)**2 * self.chi_pc + phi_s * (1-phi_s-phi_p) * (self.chi_ps + self.chi_pc - self.chi_sc) ))
	mu_c = lambda self, phi_s, phi_p: np.float64(np.log(1-phi_s-phi_p) + 1 - (1-phi_s-phi_p) - self.vc/self.vs * phi_s - self.vc/self.vp * phi_p + self.vc * (phi_s**2 * self.chi_sc + phi_p**2 * self.chi_pc + phi_s * phi_p * (self.chi_sc + self.chi_pc - self.chi_ps) ))

	delta_mu_s = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(self.mu_s(phi_s1, phi_p1) - self.mu_s(phi_s2, phi_p2))
	delta_mu_p = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(self.mu_p(phi_s1, phi_p1) - self.mu_p(phi_s2, phi_p2))
	delta_mu_c = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(self.mu_c(phi_s1, phi_p1) - self.mu_c(phi_s2, phi_p2))

	d_delta_mu_s_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-1 + 1/phi_s1 - self.chi_pc * phi_p1 * self.vs + self.chi_ps * phi_p1 * self.vs + self.vs/self.vc + self.chi_sc * ( phi_p1 - 2 * (1 - phi_p1 - phi_s1)) * self.vs)
	d_delta_mu_s_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( 1 - 1/phi_s2 + self.chi_pc * phi_p2 * self.vs - self.chi_ps * phi_p2 * self.vs - self.vs/self.vc + self.chi_sc * (-phi_p2 + 2 * (1 - phi_p2 - phi_s2)) * self.vs)
	d_delta_mu_s_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-self.chi_pc * phi_s1 * self.vs + self.vs/self.vc - self.vs/self.vp + self.chi_ps * ( 2 * phi_p1 + phi_s1) * self.vs + self.chi_sc * (-2 * (1 - phi_p1 - phi_s1) + phi_s1) * self.vs)
	d_delta_mu_s_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( self.chi_pc * phi_s2 * self.vs - self.vs/self.vc + self.vs/self.vp + self.chi_ps * (-2 * phi_p2 - phi_s2) * self.vs + self.chi_sc * ( 2 * (1 - phi_p2 - phi_s2) - phi_s2) * self.vs)

	d_delta_mu_p_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( self.vp/self.vc + self.chi_pc * (1 - phi_p1 - 2 * phi_s1) * self.vp + self.chi_ps * ( 1 - phi_p1) * self.vp + self.chi_sc * (2 * phi_s1 + phi_p1 - 1) * self.vp - self.vp/self.vs)
	d_delta_mu_p_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-self.vp/self.vc + self.chi_sc * (1 - phi_p2 - 2 * phi_s2) * self.vp + self.chi_ps * (-1 + phi_p2) * self.vp + self.chi_pc * (2 * phi_s2 + phi_p2 - 1) * self.vp + self.vp/self.vs)
	d_delta_mu_p_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-1 + 1/phi_p1 - self.chi_pc * phi_s1 * self.vp - self.chi_ps * phi_s1 * self.vp + self.chi_sc * phi_s1 * self.vp + self.vp/self.vc)
	d_delta_mu_p_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( 1 - 1/phi_p2 + self.chi_pc * phi_s2 * self.vp + self.chi_ps * phi_s2 * self.vp - self.chi_sc * phi_s2 * self.vp - self.vp/self.vc)

	d_delta_mu_c_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( 1 - 1/( 1 - phi_p1 - phi_s1) + self.chi_pc * phi_p1 * self.vc - self.chi_ps * phi_p1 * self.vc + self.chi_sc * ( phi_p1 + 2 * phi_s1) * self.vc - self.vc/self.vs)
	d_delta_mu_c_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-1 + 1/( 1 - phi_p2 - phi_s2) - self.chi_pc * phi_p2 * self.vc + self.chi_ps * phi_p2 * self.vc + self.chi_sc * (-phi_p2 - 2 * phi_s2) * self.vc + self.vc/self.vs)
	d_delta_mu_c_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64( 1 - 1/( 1 - phi_p1 - phi_s1) - self.chi_ps * phi_s1 * self.vc + self.chi_sc * phi_s1 * self.vc + self.chi_pc * ( 2 * phi_p1 + phi_s1) * self.vc - self.vc/self.vp)
	d_delta_mu_c_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: np.float64(-1 + 1/( 1 - phi_p2 - phi_s2) + self.chi_ps * phi_s2 * self.vc - self.chi_sc * phi_s2 * self.vc + self.chi_pc * (-2 * phi_p2 - phi_s2) * self.vc + self.vc/self.vp)

	def find_solution_in_nbrhd(self, crit_point, eps=0.01):

		delta_pp2 = [] 
		phi_s1 = [crit_point[0]]
		phi_s2 = [crit_point[0]]
		phi_p1 = [crit_point[1]]
		phi_p2 = [crit_point[1]]

		good_root = False
		pert = [-eps, 0, eps]
		pert_p2 = [-eps, eps]
		for dp_s1 in pert:
			for dp_p1 in pert:
				for dp_s2 in pert:
					for dp_p2 in pert_p2:
						def dmu_init(phi_):
							eq1 = B.delta_mu_s(phi_[0], phi_[1], phi_[2], crit_point[1]+dp_p2)/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							eq2 = B.delta_mu_p(phi_[0], phi_[1], phi_[2], crit_point[1]+dp_p2)/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							eq3 = B.delta_mu_c(phi_[0], phi_[1], phi_[2], crit_point[1]+dp_p2)/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							return [eq1, eq2, eq3]

						print(f"Guess provided: phi1 = ({crit_point[0]+dp_s1, crit_point[1]+dp_p1}), phi2 = ({crit_point[0]+dp_s2, crit_point[1]+dp_p2})...", flush=True)
						root = fsolve(dmu_init, [crit_point[0]+dp_s1, crit_point[1]+dp_p1, crit_point[0]+dp_s2], xtol=1e-16)
						p1   = np.array([root[0], root[1]])
						p2   = np.array([root[2], crit_point[1]+dp_p2])

						if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0 or root[2] > 1 or root[2] < 0:
							print("Breaking out...")
							continue
						#####
						if (np.abs(dmu_init(root))>1e-12).any():
							print(f"Bad root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], crit_point[1]+dp_p2})...", flush=True)
							continue
						elif np.linalg.norm(p1-p2) < 1e-3: 
							continue
						else:
							print(f"Found root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], crit_point[1]+dp_p2})...", flush=True)
							good_root = True
							delta_pp2.append(dp_p2)
							phi_s1.append(root[0])
							phi_p1.append(root[1])
							phi_s2.append(root[2])
							phi_p2.append(crit_point[1]+dp_p2)
							break

					# at this line, I am out of dp_p2 and in dp_s2
					if good_root:
						break
				# at this line, I am out of dp_s2 and in dp_p1
				if good_root:
					break
			# at this line, I am out of dp_p1 and in dp_s1
			if good_root:
				break

		if not good_root:
			print(f"No neighboring solution found...")
			exit()

		print("I have found the initial root outside of the critical point. Time to go beyond.", flush=True)
		print(f"phi_s1 = {phi_s1}", flush=True)
		print(f"phi_p1 = {phi_p1}", flush=True)
		print(f"phi_s2 = {phi_s2}", flush=True)
		print(f"phi_p2 = {phi_p2}", flush=True)

		return phi_s1, phi_s2, phi_p1, phi_p2, delta_pp2

	def calc_perts(self, phi1, phi2, delta_phip2):

		a1 = np.float64(self.d_delta_mu_s_dps1(phi1[0], phi1[1], phi2[0], phi2[1]))
		a2 = np.float64(self.d_delta_mu_s_dpp1(phi1[0], phi1[1], phi2[0], phi2[1]))
		a3 = np.float64(self.d_delta_mu_s_dps2(phi1[0], phi1[1], phi2[0], phi2[1]))
		a4 = np.float64(-delta_phip2 * self.d_delta_mu_s_dpp2(phi1[0], phi1[1], phi2[0], phi2[1]))

		b1 = np.float64(self.d_delta_mu_p_dps1(phi1[0], phi1[1], phi2[0], phi2[1]))
		b2 = np.float64(self.d_delta_mu_p_dpp1(phi1[0], phi1[1], phi2[0], phi2[1]))
		b3 = np.float64(self.d_delta_mu_p_dps2(phi1[0], phi1[1], phi2[0], phi2[1]))
		b4 = np.float64(-delta_phip2 * self.d_delta_mu_p_dpp2(phi1[0], phi1[1], phi2[0], phi2[1]))

		c1 = np.float64(self.d_delta_mu_c_dps1(phi1[0], phi1[1], phi2[0], phi2[1]))
		c2 = np.float64(self.d_delta_mu_c_dpp1(phi1[0], phi1[1], phi2[0], phi2[1]))
		c3 = np.float64(self.d_delta_mu_c_dps2(phi1[0], phi1[1], phi2[0], phi2[1]))
		c4 = np.float64(-delta_phip2 * self.d_delta_mu_c_dpp2(phi1[0], phi1[1], phi2[0], phi2[1]))

		D  = np.float64(np.linalg.det(np.array([[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]], dtype=np.float64)))
		Dx = np.float64(np.linalg.det(np.array([[a4,a2,a3],[b4,b2,b3],[c4,c2,c3]], dtype=np.float64)))
		Dy = np.float64(np.linalg.det(np.array([[a1,a4,a3],[b1,b4,b3],[c1,c4,c3]], dtype=np.float64)))
		Dz = np.float64(np.linalg.det(np.array([[a1,a2,a4],[b1,b2,b4],[c1,c2,c4]], dtype=np.float64)))

		# print(f"a1 = {a1}, a2 = {a2}, a3 = {a3}, a4 = {a4}")
		# print(f"D = {D}, Dx = {Dx}, Dy = {Dy}, Dz = {Dz}")
		# print(f"dx = {Dx/D}, dy = {Dy/D}, dz = {Dz/D}")

		return Dx/D, Dy/D, Dz/D

	def run_cycles(self, ncycles, crit_point, continuation_b):
		if continuation_b:
			pass
		else:
			phi_s1, phi_s2, phi_p1, phi_p2, delta_pp2 = self.find_solution_in_nbrhd_pert_p2(crit_point)

		print(f"delta_pp2 = {delta_pp2}")
		condition = (phi_s1[-1] < 1e-3 or phi_p1[-1] < 1e-3 or 1-phi_s1[-1]-phi_p1[-1] < 1e-3)

		# for i in range(ncycles):
		iterr  = 0 
		max_it = 1e+4

		while not condition:
			iterr += 1
			if iterr > max_it:
				break

			print(f"@ i = {iterr}/{max_it}...", flush=True)

			phi1 = [phi_s1[-1], phi_p1[-1]]
			phi2 = [phi_s2[-1], phi_p2[-1]]
			delta_ps1, delta_pp1, delta_ps2 = self.calcs_perts(phi1, phi2, delta_pp2[-1]) # delta_finder_pert_p2(phi1, phi2, delta_pp2[-1])
			print(f"delta_ps1 = {delta_ps1}, delta_pp1 = {delta_pp1}, delta_ps2 = {delta_ps2}, delta_pp2 = {delta_pp2[-1]}", flush=True)

			def dmu(phi_):
				eq1 = self..delta_mu_s(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2[-1])
				eq2 = self.delta_mu_p(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2[-1])
				eq3 = self.delta_mu_c(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2[-1])
				return [eq1, eq2, eq3]

			print(f"Guess provided: phi1 = ({phi_s1[-1]+delta_ps1, phi_p1[-1]+delta_pp1}), phi2 = {phi_s2[-1]+delta_ps2, phi_p2[-1]+delta_pp2[-1]}", flush=True)
			root = fsolve(dmu, [phi_s1[-1]+delta_ps1, phi_p1[-1]+delta_pp1, phi_s2[-1]+delta_ps2], xtol=1e-16)

			if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0 or root[2] > 1 or root[2] < 0:
				print("Breaking out...")
				break

			p1 = np.array([root[0], root[1]])
			p2 = np.array([root[2], phi_p2[-1]+delta_pp2[-1]])

			p1_ = np.array([phi_s1[-1], phi_p1[-1]])
			p2_ = np.array([phi_s2[-1], phi_p2[-1]])

			d1 = p1 - p1_ 
			d2 = p2 - p2_ 

			p1__ = np.array([phi_s1[-2], phi_p1[-2]])
			p2__ = np.array([phi_s2[-2], phi_p2[-2]])

			d1_ = p1_ - p1__
			d2_ = p2_ - p2__


			if np.linalg.norm(root[0] - root[2]) > 1e-4:
				print(f"p1 = {p1}, p2 = {p2}")
				print(f"p1_ = {p1_}, p2_ = {p2_}")
				print(f"p1__ = {p1__}, p2__ = {p2__}")
				print ("PROBLEM!")

			if (np.abs(dmu(root))>1e-12).any():
				print(f"Bad root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], phi_p2[-1]+delta_pp2[-1]})...", flush=True)
				delta_pp2.append(delta_pp2[-1]/2)
				continue
			elif np.linalg.norm(p1-p2) < 1e-3: 
				print(f"Too close: phi1 = ({root[0], root[1]}), phi2 = ({root[2], phi_p2[-1]+delta_pp2[-1]})...", flush=True)
				delta_pp2.append(delta_pp2[-1]/2)
				continue
			elif np.dot(d1/np.linalg.norm(d1), d1_/np.linalg.norm(d1_)) <= 0 or np.dot(d2/np.linalg.norm(d2), d2_/np.linalg.norm(d2_)) <= 0:
				print(f"Making an about turn: phi1 = ({root[0], root[1]}), phi2 = ({root[2], phi_p2[-1]+delta_pp2[-1]})...", flush=True)
				delta_pp2.append(delta_pp2[-1]/2)
				continue


			print(f"Found root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], phi_p2[-1]+delta_pp2[-1]})...", flush=True)

			phi_s1.append(root[0])
			phi_p1.append(root[1])
			phi_s2.append(root[2])
			phi_p2.append(phi_p2[-1]+delta_pp2[-1])
			delta_pp2.append(delta_pp2[-1])
			condition = (phi_s1[-1] < 1e-3 or phi_p1[-1] < 1e-3 or 1-phi_s1[-1]-phi_p1[-1] < 1e-3)

		phi_s1 = np.array(phi_s1)
		phi_p1 = np.array(phi_p1)
		phi_s2 = np.array(phi_s2)
		phi_p2 = np.array(phi_p2)

		return phi_s1, phi_s2, phi_p1, phi_p2

# end of class Mu_PS

class Mu_PC:
	def __init__(self, inputs):
		self.chi_sc = inputs["chi_sc"]
		self.chi_ps = inputs["chi_ps"]
		self.chi_pc = inputs["chi_pc"]
		self.vs     = inputs["vs"]
		self.vc     = inputs["vc"]
		self.vp     = inputs["vp"]
		return

	# calculate chemical potentials
	mu_s = lambda self, phi_p, phi_c: np.float64(np.log(1-phi_p-phi_c) + 1 - (1 - phi_p - phi_c) - self.vs/self.vp * phi_p               - self.vs/self.vc * phi_c + self.vs * (phi_p**2               * self.chi_ps + phi_c**2 * self.chi_sc + phi_p * (1 - phi_p - phi_c) * (self.chi_ps + self.chi_sc - self.chi_pc)))
	mu_p = lambda self, phi_p, phi_c: np.float64(np.log(phi_p)         + 1 - phi_p               - self.vp/self.vs * (1 - phi_p - phi_c) - self.vp/self.vc * phi_c + self.vp * ((1 - phi_p - phi_c)**2 * self.chi_sc + phi_p**2 * self.chi_pc + (1 - phi_c - phi_p) * phi_p * (self.chi_ps + self.chi_pc - self.chi_sc)))
	mu_c = lambda self, phi_p, phi_c: np.float64(np.log(phi_c)         + 1 - phi_c               - self.vc/self.vp * phi_p               - self.vc/self.vp * phi_s + self.vc * ((1 - phi_p - phi_c)**2 * self.chi_sc + phi_p**2 * self.chi_pc + (1 - phi_c - phi_p) * phi_p * (self.chi_sc + self.chi_pc - self.chi_ps)))

	delta_mu_s = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: np.float64(self.mu_s(phi_p1, phi_c1) - self.mu_s(phi_p2, phi_c2))
	delta_mu_p = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: np.float64(self.mu_p(phi_p1, phi_c1) - self.mu_p(phi_p2, phi_c2))
	delta_mu_c = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: np.float64(self.mu_c(phi_p1, phi_c1) - self.mu_c(phi_p2, phi_c2))

	d_delta_mu_s_dpp1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: 1  - 1/(1 - phi_p1 - phi_c1) - self.vs/self.vp + self.chi_sc * ((1 - phi_c1 - phi_p1) - phi_p1) * self.vs + self.chi_ps * ( (1 - phi_c1 - phi_p1) + phi_p1) * self.vs + self.chi_pc * (phi_p1 + (-1 + phi_c1 + phi_p1)) * self.vs
	d_delta_mu_s_dpp2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: -1 + 1/(1 - phi_p2 - phi_c2) + self.vs/self.vp + self.chi_pc * ((1 - phi_c2 - phi_p2) - phi_p2) * self.vs + self.chi_ps * (-phi_p2 + (-1 +phi_c2 + phi_p2)) * self.vs + self.chi_sc * (phi_p2 + (-1 + phi_c2 + phi_p2)) * self.vs
	d_delta_mu_s_dpc1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: 1  - 1/(1 - phi_c1 - phi_p1) + self.chi_pc * phi_p1 * self.vs - self.chi_ps * phi_p1 * self.vs - self.vs/self.vc + self.chi_sc * ( 2 * phi_c1 - phi_p1) * self.vs
	d_delta_mu_s_dpc2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: -1 + 1/(1 - phi_c2 - phi_p2) - self.chi_pc * phi_p2 * self.vs + self.chi_ps * phi_p2 * self.vs + self.vs/self.vc + self.chi_sc * (-2 * phi_cs + phi_p2) * self.vs

	d_delta_mu_p_dpp1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: -1 + 1/phi_p1 - self.chi_pc * phi_c1 * self.vp + self.chi_sc * phi_c1 * self.vp + self.chi_ps * (-phi_c1 - 2 * (1 - phi_c1 - phi_p1)) * self.vp + self.vp/self.vs
	d_delta_mu_p_dpp2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2:  1 - 1/phi_p1 + self.chi_pc * phi_c2 * self.vp - self.chi_sc * phi_c2 * self.vp + self.chi_ps * ( phi_c2 + 2 * (1 - phi_c2 - phi_p2)) * self.vp - self.vp/self.vs
	d_delta_mu_p_dpc1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: self.chi_pc * (1 - phi_p1) * self.vp + self.chi_ps * (-1 + phi_p1) * self.vp             + self.chi_sc * (-1 + 2 * phi_c1 + phi_p1) * self.vp + self.vp * (-1/self.vc + 1/self.vs)
	d_delta_mu_p_dpc2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: self.chi_ps * (1 - phi_p2) * self.vp + self.chi_sc * (1 - 2 * phi_c2 - phi_p2) * self.vp + self.chi_pc * (-1 + phi_p2) * self.vp              + self.vp * ( 1/self.vc - 1/self.vs)

	d_delta_mu_c_dpp1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: self.chi_sc * (-1 + phi_c1) * self.vc + self.chi_ps * (-1 + phi_c1 + 2 * phi_p1) * self.vc + self.chi_pc * (1 - phi_c1) * self.vc            - self.vc/self.vp + self.vc/self.vs
	d_delta_mu_c_dpp2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: self.chi_sc * ( 1 - phi_c2) * self.vc + self.chi_pc * (-1 + phi_c2) * self.vc              + self.chi_ps * (1 - phi_c2 - 2*phi_p2) * self.vc + self.vc/self.vp - self.vc/self.vp
	d_delta_mu_c_dpc1 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2: -1 + 1/phi_c1 - self.chi_pc * phi_p1 * self.vc + self.chi_ps * phi_p1 * self.vc + self.chi_sc * ( 2 * phi_c1 + (-2 + phi_p1)) * self.vc + self.vc/self.vs
	d_delta_mu_c_dpc2 = lambda self, phi_p1, phi_c1, phi_p2, phi_c2:  1 - 1/phi_c2 + self.chi_pc * phi_p2 * self.vc - self.chi_ps * phi_p2 * self.vc + self.chi_sc * (-2 * phi_c2 + ( 2 - phi_p2)) * self.vc - self.vc/self.vs

	def find_solution_in_nbrhd(self, crit_point, eps=0.01):

		delta_pp2 = [] 
		phi_p1 = [crit_point[0]]
		phi_c1 = [crit_point[0]]
		phi_p2 = [crit_point[1]]
		phi_c2 = [crit_point[1]]

		good_root = False
		pert = [-eps, 0, eps]
		pert_p2 = [-eps, eps]
		for dp_c1 in pert:
			for dp_p1 in pert:
				for dp_c2 in pert:
					for dp_p2 in pert_p2:
						def dmu_init(phi_):
							eq1 = B.delta_mu_s(phi_[0], phi_[1], crit_point[0]+dp_p2, phi_[2])/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							eq2 = B.delta_mu_p(phi_[0], phi_[1], crit_point[0]+dp_p2, phi_[2])/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							eq3 = B.delta_mu_c(phi_[0], phi_[1], crit_point[0]+dp_p2, phi_[2])/np.linalg.norm(np.array([phi_[0], phi_[1]]) - np.array([phi_[2], crit_point[1]+dp_p2]))
							return [eq1, eq2, eq3]

						print(f"Guess provided: phi1 = ({crit_point[0]+dp_s1, crit_point[1]+dp_p1}), phi2 = ({crit_point[0]+dp_s2, crit_point[1]+dp_p2})...", flush=True)
						root = fsolve(dmu_init, [crit_point[0]+dp_s1, crit_point[1]+dp_p1, crit_point[0]+dp_s2], xtol=1e-16)
						p1   = np.array([root[0], root[1]])
						p2   = np.array([root[2], crit_point[1]+dp_p2])

						if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0 or root[2] > 1 or root[2] < 0:
							print("Breaking out...")
							continue
						#####
						if (np.abs(dmu_init(root))>1e-12).any():
							print(f"Bad root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], crit_point[1]+dp_p2})...", flush=True)
							continue
						elif np.linalg.norm(p1-p2) < 1e-3: 
							continue
						else:
							print(f"Found root: phi1 = ({root[0], root[1]}), phi2 = ({root[2], crit_point[1]+dp_p2})...", flush=True)
							good_root = True
							delta_pp2.append(dp_p2)
							phi_s1.append(root[0])
							phi_p1.append(root[1])
							phi_s2.append(root[2])
							phi_p2.append(crit_point[1]+dp_p2)
							break

					# at this line, I am out of dp_p2 and in dp_s2
					if good_root:
						break
				# at this line, I am out of dp_s2 and in dp_p1
				if good_root:
					break
			# at this line, I am out of dp_p1 and in dp_s1
			if good_root:
				break

		if not good_root:
			print(f"No neighboring solution found...")
			exit()

		print("I have found the initial root outside of the critical point. Time to go beyond.", flush=True)
		print(f"phi_s1 = {phi_s1}", flush=True)
		print(f"phi_p1 = {phi_p1}", flush=True)
		print(f"phi_s2 = {phi_s2}", flush=True)
		print(f"phi_p2 = {phi_p2}", flush=True)

		return phi_s1, phi_s2, phi_p1, phi_p2, delta_pp2


