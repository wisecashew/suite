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
	mu_s = lambda self, phi_s, phi_p: np.log(phi_s)         + 1 - phi_s - self.vs/self.vp * phi_p - self.vs/self.vc * (1-phi_s-phi_p) + self.vs * (phi_p**2 * self.chi_ps + (1-phi_s-phi_p)**2 * self.chi_sc + phi_p * (1-phi_s-phi_p) * (self.chi_ps + self.chi_sc - self.chi_pc) )
	mu_p = lambda self, phi_s, phi_p: np.log(phi_p)         + 1 - phi_p - self.vp/self.vs * phi_s - self.vp/self.vc * (1-phi_s-phi_p) + self.vp * (phi_s**2 * self.chi_ps + (1-phi_s-phi_p)**2 * self.chi_pc + phi_s * (1-phi_s-phi_p) * (self.chi_ps + self.chi_pc - self.chi_sc) )
	mu_c = lambda self, phi_s, phi_p: np.log(1-phi_s-phi_p) + 1 - (1-phi_s-phi_p) - self.vc/self.vs * phi_s - self.vc/self.vp * phi_b + self.vc * (phi_s**2 * self.chi_sc + phi_p**2 * self.chi_pc + phi_s * phi_p * (self.chi_sc + self.chi_pc - self.chi_ps) )

	delta_mu_s = lambda self, phi_s1, phi_b1, phi_a2, phi_b2: self.mu_a(phi_s1, phi_p1) - self.mu_a(phi_s2, phi_p2) 
	delta_mu_p = lambda self, phi_s1, phi_b1, phi_a2, phi_b2: self.mu_b(phi_s1, phi_p1) - self.mu_b(phi_s2, phi_p2)
	delta_mu_c = lambda self, phi_s1, phi_b1, phi_a2, phi_b2: self.mu_c(phi_s1, phi_p1) - self.mu_c(phi_s2, phi_p2)

	d_delta_mu_s_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: -1 + 1/phi_s1 + (-2*self.chi_sc*(1-phi_s1-phi_p1) - (-self.chi_pc + self.chi_ps + self.chi_sc) * phi_p1) * self.vs + self.vs/self.vc
	d_delta_mu_s_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2:  1 - 1/phi_s2 - (-2*self.chi_sc*(1-phi_s2-phi_p2) - (-self.chi_pc + self.chi_ps + self.chi_sc) * phi_p1) * self.vs - self.vs/self.vc
	d_delta_mu_s_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: (self.chi_ps - self.chi_sc - self.chi_ps * phi_s1 + self.chi_pc * (-1 + phi_s1 + 2*phi_b1) + 1/self.vc - 1/self.vp) * self.vs
	d_delta_mu_s_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: (self.chi_sc + self.chi_ps*(-1 + phi_s2) - self.chi_sc*phi_s2 - self.chi_pc * (-1 + phi_s2 + 2*phi_p2) - 1/self.vc + 1/self.vp) * vs

	d_delta_mu_p_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: self.vp*(self.chi_ps - self.chi_sc + 2*self.chi_sc*phi_a1 + self.chi_pc * (phi_b1 - 1) - self.chi_ps * phi_b1 + self.chi_sc * phi_b1 + 1/vc - 1/vs)
	d_delta_mu_p_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: self.vp*(self.chi_pc + self.chi_sc - 2*self.chi_sc*phi_a2 + self.chi_ps * (phi_b2 - 1) - self.chi_pc * phi_b2 - self.chi_sc * phi_b2 - 1/vc + 1/vs)
	d_delta_mu_p_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: -1 + 1/phi_p1 + self.chi_pc * (-2 + phi_s1)*self.vp - self.chi_ps * phi_s1 * self.vp + self.chi_sc * phi_s1 * self.vp + 2 * self.chi_pc * phi_p1 * self.vp + self.vp/self.vc
	d_delta_mu_p_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2:  1 - 1/phi_p2 - self.chi_pc * (-2 + phi_s2)*self.vp + self.chi_ps * phi_s2 * self.vp - self.chi_sc * phi_s2 * self.vp - 2 * self.chi_pc * phi_p2 * self.vp - self.vp/self.vc

	d_delta_mu_c_dps1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2:  1 + 1/(-1 + phi_s1 + phi_p1) + 2 * self.chi_sc * phi_s1 * self.vc + (self.chi_pc - self.chi_ps + self.chi_sc) * phi_p1 * self.vc - self.vc/self.vs
	d_delta_mu_c_dps2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: -1 - 1/(-1 + phi_s2 + phi_p2) - (2 * self.chi_sc * phi_s2 + (self.chi_pc - self.chi_ps + self.chi_sc) * phi_p2) * self.vc + self.vc/self.vs
	d_delta_mu_c_dpp1 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2:  1 + 1/(-1 + phi_s1 + phi_p1) + (self.chi_pc - self.chi_ps + self.chi_sc) * phi_s1 * self.vc + 2 * self.chi_pc * phi_p1 * self.vc - self.vc/self.vp
	d_delta_mu_c_dpp2 = lambda self, phi_s1, phi_p1, phi_s2, phi_p2: -1 - 1/(-1 + phi_s2 + phi_p2) - ((self.chi_pc - self.chi_ps + self.chi_sc) * phi_s2 + 2 * self.chi_pc * phi_p2) * vc + self.vc/self.vp

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

	def D_det_calcs(self, phi1, phi2, delta_phip2):
		a1 = self.d_delta_mu_s_dps1(phi1[0], phi1[1], phi2[0], phi2[1])
		a2 = self.d_delta_mu_s_dpp1(phi1[0], phi1[1], phi2[0], phi2[1])
		a3 = self.d_delta_mu_s_dps2(phi1[0], phi1[1], phi2[0], phi2[1])
		a4 = -delta_phip2 * self.d_delta_mu_s_dpp2(phi1[0], phi1[1], phi2[0], phi2[1])

		b1 = self.d_delta_mu_p_dps1(phi1[0], phi1[1], phi2[0], phi2[1])
		b2 = self.d_delta_mu_p_dpp1(phi1[0], phi1[1], phi2[0], phi2[1])
		b3 = self.d_delta_mu_p_dps2(phi1[0], phi1[1], phi2[0], phi2[1])
		b4 = -delta_phip2 * self.d_delta_mu_p_dpp2(phi1[0], phi1[1], phi2[0], phi2[1])

		c1 = self.d_delta_mu_c_dps1(phi1[0], phi1[1], phi2[0], phi2[1])
		c2 = self.d_delta_mu_c_dps1(phi1[0], phi1[1], phi2[0], phi2[1])
		c3 = self.d_delta_mu_c_dps1(phi1[0], phi1[1], phi2[0], phi2[1])
		c4 = -delta_phip2 * self.d_delta_mu_c_dps1(phi1[0], phi1[1], phi2[0], phi2[1])

		D  = np.array([[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]])
		Dx = np.array([[a4,a2,a3],[b4,b2,b3],[c4,c2,c3]])
		Dy = np.array([[a1,a4,a3],[b1,b4,b3],[c1,c4,c3]])
		Dz = np.array([[a1,a2,a4],[b1,b2,b4],[c1,c2,c4]])

		return np.linalg.det(D), np.linalg.det(Dx), np.linalg.det(Dy), np.linalg.det(Dz)

	def delta_finder(self, phi1, phi2, delta_phip2):

		D, Dx, Dy, Dz  = self.D_det_calcs (self, phi1, phi2, delta_phip2)

		delta_phis1 = Dx/D
		delta_phip1 = Dy/D
		delta_phis2 = Dz/D

		return [delta_phis1, delta_phip1, delta_phis2]

	######################################################
	# start performing sweeps
	# first sweep

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

	# define the center and crit point around which you want to find the binodal
	center     = np.mean(B.crits, axis=0)
	crit_point = B.crits[0]
	print(f"crit_point = {crit_point}", flush=True)

	#######
	delta_pp2 = 0.01 
	phi_s1 = [crit_point[0]]
	phi_s2 = [crit_point[0]]
	phi_p1 = [crit_point[1]]
	phi_p2 = [crit_point[1]]

	ncycles = int(crit_point[1]/delta_pb2) + 1

	for i in range(ncycles):
		print(f"@ i = {i}/{ncycles}...", flush=True)

		if i > 0:
			phi1 = [phi_s1[-1], phi_p1[-1]]
			phi2 = [phi_s2[-1], phi_p2[-1]]
			deltas = B.delta_finder(phi1, phi2, delta_pp2)
			delta_ps1 = deltas[0]
			delta_pp1 = deltas[1]
			delta_ps2 = deltas[2]
		else:
			delta_ps1 = -0.01
			delta_pp1 = -0.01
			delta_ps2 = 0.01

		def dmu(phi_):
			eq1 = self.delta_mu_s(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2)
			eq2 = self.delta_mu_p(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2)
			eq3 = self.delta_mu_c(phi_[0], phi_[1], phi_[2], phi_p2[-1]+delta_pp2)
			return [eq1, eq2, eq3]

		root = fsolve(dmu, [phi_s1[-1]+delta_ps1, phi_p1[-1]+delta_pp1, phi_s2[-1]+delta_ps2], xtol=1e-9)

		if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0 or root[2] > 1 or root[2] < 0:
			print("Breaking out...")
			break

		if (abs(dmu(root))>1e-9).any()
			continue

		phi_s1.append(root[0])
		phi_p1.append(root[1])
		phi_s2.append(root[2])
		phi_p2.append(phi_p2[-1]+delta_pp2)

	phi_s1 = np.array(phi_s1)
	phi_p1 = np.array(phi_p1)
	phi_s2 = np.array(phi_s2)
	phi_p2 = np.array(phi_p2)

	ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='darkred',   s=1)
	ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='slategray', s=1)

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

	stop = time.time()
	print (f"Time taken to scan the ternary space has been {stop-start} seconds.", flush=True)
