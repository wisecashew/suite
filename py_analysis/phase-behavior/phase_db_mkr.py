import numpy as np 
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
import argparse 

#####################################################
# create a parser to take in the inputs

parser = argparse.ArgumentParser(description="Creates a database of points on the T-phi diagram and their thermodynamic state per the parameters you enter (set kb = 1).")
parser.add_argument("--EMM-A",    dest='emma', type=float, action='store', help="Value of aligned monomer-monomer interaction (required).")
parser.add_argument("--EMM-N",    dest='emmn', type=float, action='store', help="Value of misaligned monomer-monomer interaction (required).")
parser.add_argument("--EMS-A",    dest='emsa', type=float, action='store', help="Value of aligned monomer-solvent interaction (required).")
parser.add_argument("--EMS-N",    dest='emsn', type=float, action='store', help="Value of misaligned monomer-solvent interaction (required).")
parser.add_argument("--ESS-A",    dest='essa', type=float, action='store', help="Value of aligned solvent-solvent interaction (required).")
parser.add_argument("--ESS-N",    dest='essn', type=float, action='store', help="Value of misaligned solvent-solvent interaction (required).")
parser.add_argument("--PV",       dest='pv',   type=float, action='store', help="Value of operational volume (required, default: 0).", default=0)
parser.add_argument("--PW",       dest='pw',   type=float, action='store', help="Value of microcanonical probability of being aligned (required, default: 0).", default=0)
parser.add_argument("-N",         dest='n',    type=int,   action='store', help="Enter degree of polymerization of your polymer (required, default: 32).", default=32)
parser.add_argument("--phi-mesh", dest='pm',   type=int,   action='store', help="Enter mesh size of phi points from 0 to 1 (default: 1000).", default=1000)
parser.add_argument("--T-mesh",   dest='tm',   type=int,   action='store', help="Enter mesh size of T points from 0.001 to 100 (default: 1000).", default=1000)
parser.add_argument("--dumpfile", dest='df',   type=str,   action='store', help="Enter name of file to dump out (phi, T, state) information.")
parser.add_argument("--to-plot",  dest='plot', action='store_true', help="Enter this option if you want to plot the spinodal (default: False)", default=False)
parser.add_argument("--img-name", dest='img',  type=str,   action='store', help="Enter name of image file to be created (dont put the .png in the input) (default: spinodal).", default="spinodal")
args = parser.parse_args()

# end of argument parser
#####################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################
# utility: in-line auxiliary functions for chi
zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float128)
zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float128)
zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128) + (1-pw)*np.exp (-1/T * essn, dtype=np.float128)
fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128)/zmm(emma, emmn, pw, T)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128)/zms(emsa, emsn, pw, T)
fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128)/zss(essa, essn, pw, T)

#####################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################
# create a Phase class that holds all the information to run the calculations
# this helps to keep the code compact and organized

class Phase:

	# utility: initializes the phase object
	# dependencies: an dictionary "inputs" that has all the information about the energetic parameters
	def __init__(self, inputs):
		self.emm_a = inputs["emm_a"]
		self.emm_n = inputs["emm_n"]
		self.ems_a = inputs["ems_a"]
		self.ems_n = inputs["ems_n"]
		self.ess_a = inputs["ess_a"]
		self.ess_n = inputs["ess_n"]
		self.pv    = inputs["pv"]
		self.pw    = inputs["pw"]
		self.n     = inputs["n"]
		return

	# utility: function that computes chi given a list of temperatures
	# dependencies: a list of temperatures to calculate chi at 
	def chi (self, T):
		c = 24 * (self.pv * ( (fmsa (self.ems_a, self.ems_n, self.pw, T) * self.ems_a + (1 - fmsa (self.ems_a, self.ems_n, self.pw, T) ) * self.ems_n) - 1/2 * \
		( (fmma (self.emm_a, self.emm_n, self.pw, T) * self.emm_a + (1-fmma (self.emm_a, self.emm_n, self.pw, T) ) * self.emm_n) + \
		  (fssa (self.ess_a, self.ess_n, self.pw, T) * self.ess_a + (1-fssa (self.ess_a, self.ess_n, self.pw, T) ) * self.ess_n) ) )
		+ (1-self.pv) * (self.ems_n - 1/2 * (self.emm_n + self.ess_n) ) ) / T

		return np.array(c, dtype=np.float128)

	# utility: calculates the spinodal of the system
	# dependence: a list of temperatures to calculate the spinodal at
	def spinodal (self, T):
		p1 = -1/(4 * self.n * self.chi (T)) * (-1 + self.n - 2 * self.n * self.chi (T) - np.sqrt (-8 * self.n * self.chi (T) + (1 - self.n + 2 * self.n * self.chi (T) ) ** 2 ) )
		p2 = -1/(4 * self.n * self.chi (T)) * (-1 + self.n - 2 * self.n * self.chi (T) + np.sqrt (-8 * self.n * self.chi (T) + (1 - self.n + 2 * self.n * self.chi (T) ) ** 2 ) )

		return (np.array(p1, dtype=np.float128), np.array(p2, dtype=np.float128), T)

	# utility: determines if the point is homogeneous or not
	# dependence: volume fraction of polymer and a temperature
	def state(self, phi, T):
		condition = 1/(self.n*phi) + 1/(1-phi) - 2*self.chi(T)
		condition = np.where(condition>0, 1, -1)
		return condition

# end of Phase class
#####################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__=="__main__":

	start = time.time()
	system_inputs = dict()
	try:
		system_inputs["emm_a"] = args.emma
		system_inputs["emm_n"] = args.emmn
		system_inputs["ems_a"] = args.emsa
		system_inputs["ems_n"] = args.emsn
		system_inputs["ess_a"] = args.essa
		system_inputs["ess_n"] = args.essn
		system_inputs["pv"]    = args.pv
		system_inputs["pw"]    = args.pw
		system_inputs["n"]     = args.n
	except:
		print("You missed some inputs. Do `python phase_db_mkr.py -h` to check. Exiting...")
		exit()

	print("Obtained all the inputs. Making phase object...")

	# make the Phase object
	phase = Phase(system_inputs) 

	print("Made Phase object! Making grids...")

	# define the grids
	phi_grid = np.linspace(0,  1, args.pm)
	T_grid   = np.logspace(-3, 2, args.tm)

	print("Made grids! Making the mesh...")

	# make the mesh 
	phi_mesh, T_mesh = np.meshgrid(phi_grid, T_grid)

	print("Made mesh! Getting the state labels...")

	# get the state of the system
	state_labels = phase.state(phi_mesh, T_mesh)

	# plot the colormesh 
	if args.plot:
		print("Plotting the colormesh...")
		fig = plt.figure(num=0, figsize=(3,3))
		ax  = plt.axes()

		ax.pcolormesh(phi_mesh, T_mesh, state_labels, cmap='bwr')
		fig.savefig(args.img+"_mesh.png", dpi=1200, bbox_inches="tight")




	print("Got the state labels! Writing it all out...")

	# flatten the arrays out
	phi_mesh     = phi_mesh.flatten()
	T_mesh       = T_mesh.flatten()
	state_labels = state_labels.flatten()

	# start writing the stuff
	f = open(args.df, 'w')
	for i in range(len(phi_mesh)):
		f.write(f"{phi_mesh[i]}\t{T_mesh[i]}\t{state_labels[i]}\n")
	f.close()





	if args.plot:
		print("Plotting the spinodal...")
		fig = plt.figure(num=1, figsize=(3,3))
		ax  = plt.axes()
		ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')

		# mechanics...
		spinodal = phase.spinodal (T_grid)
		T1       = spinodal[2][spinodal[0] < 1]
		arm1     = spinodal[0][spinodal[0] < 1]
		T1       = T1   [arm1 > 0]
		arm1     = arm1 [arm1 > 0]
		line = ax.plot  (arm1, T1, lw=1.0, markersize=0, c="darkred", solid_capstyle="round",label="_nolabel_")
		T1       = spinodal[2][spinodal[1] < 1]
		arm1     = spinodal[1][spinodal[1] < 1]
		T1       = T1   [arm1 > 0]
		arm1     = arm1 [arm1 > 0]
		ax.plot  (arm1, T1, lw=1.0, markersize=0, color="steelblue", solid_capstyle="round", label="_nolabel_")

		fig.savefig(args.img+"_line.png", dpi=1200, bbox_inches="tight")



	end = time.time()
	print(f"Time for execution is {end-start} seconds.")
