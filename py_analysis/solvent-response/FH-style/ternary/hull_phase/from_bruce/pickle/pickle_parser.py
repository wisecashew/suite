import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
import mu
import spinodal
import mpltern
import linecache
import warnings
from scipy.optimize import fsolve

###################################################
def custom_warning_format(message, category, filename, lineno, line=None):
	line = linecache.getline(filename, lineno).strip()
	return "" # return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
###################################################

# set up the description
parser = argparse.ArgumentParser(description="Parse the newton input pickle file.")
parser.add_argument("-i", metavar="pickle file", type=str, help="Enter address of pickle file.")
parser.add_argument("-o", metavar="image file", type=str, help="Enter name of image to be generated.")
args = parser.parse_args()

if __name__=="__main__":

	fig = plt.figure(figsize=(8,8))
	ax  = fig.add_subplot(projection="ternary")
	ax.set_tlabel ("$\\phi _{S}$")
	ax.set_llabel ("$\\phi _{C}$")
	ax.set_rlabel ("$\\phi _{P}$")
	ax.set_tlim(0, 1)
	ax.set_llim(0, 1)
	ax.set_rlim(0, 1)
	positions = ['tick1', 'tick2']
	for position in positions:
		ax.taxis.set_ticks_position(position)
		ax.laxis.set_ticks_position(position)
		ax.raxis.set_ticks_position(position)

	ax.taxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.raxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.laxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	
	f = open(args.i, 'rb')
	NEWTON = pickle.load(f)
	f.close()
	L_newton = NEWTON.shape[0]

	# get the inputs
	x_input = NEWTON[:, :8]
	y_class_true = NEWTON[:, 8]
	y_class_pred = NEWTON[:, 9]
	# print(y_class_pred)
	y_comp_true  = NEWTON[:, 10:19]
	y_comp_pred  = NEWTON[:, 19:28]
	phase_idx    = NEWTON[:, -1]

	# set up the input dictionary 
	inputs = dict()
	inputs["vs"]     = x_input[0,0]
	inputs["vc"]     = x_input[0,1]
	inputs["vp"]     = x_input[0,2]
	inputs["chi_sc"] = x_input[0,3]
	inputs["chi_ps"] = x_input[0,4]
	inputs["chi_pc"] = x_input[0,5]
	inputs["phi_s"]  = x_input[0,6]
	inputs["phi_p"]  = x_input[0,7]
	S = spinodal.Spinodal(inputs)
	M = mu.sym_mu_ps(inputs, S)

	# get triangular points
	def triangle_mu(P):
		eq1 = M.delta_mu_s(P[0], P[1], P[4], P[5])
		eq2 = M.delta_mu_p(P[0], P[1], P[4], P[5])
		eq3 = M.delta_mu_c(P[0], P[1], P[4], P[5])
		eq4 = M.delta_mu_s(P[2], P[3], P[4], P[5])
		eq5 = M.delta_mu_p(P[2], P[3], P[4], P[5])
		eq6 = M.delta_mu_c(P[2], P[3], P[4], P[5])
		return [eq1, eq2, eq3, eq4, eq5, eq6]

	net_line = 0
	bad_line = 0
	net_triangle = 0
	bad_triangle = 0
	for i in range(L_newton):
		if y_class_pred[i] == 0:
			pass
		elif y_class_pred[i] == 1:
			net_line += 1
			def line_mu(P, phi_s):
				eq1 = M.delta_mu_s(phi_s, P[0], P[1], P[2])
				eq2 = M.delta_mu_p(phi_s, P[0], P[1], P[2])
				eq3 = M.delta_mu_c(phi_s, P[0], P[1], P[2])
				return [eq1, eq2, eq3]
			ps1 = y_comp_pred[i][0]
			pp1 = y_comp_pred[i][1]
			ps2 = y_comp_pred[i][2]
			pp2 = y_comp_pred[i][3]
			root = fsolve(line_mu, [pp1, ps2, pp2], args=(ps1))
			error = line_mu(root, ps1)
			# print(f"line error = {line_mu(root, ps1)}", flush=True)
			p1 = np.array([ps1, root[0]])
			p2 = np.array([root[1], root[2]])
			if (np.abs(error) < 1e-6).all(): 
				ax.plot([ps1, root[1]], [1-ps1-root[0],1-root[1]-root[2]], [root[0], root[2]], c="coral", markersize=0.5)
			def line_mu(P, phi_p):
				eq1 = M.delta_mu_s(P[0], phi_p, P[1], P[2])
				eq2 = M.delta_mu_p(P[0], phi_p, P[1], P[2])
				eq3 = M.delta_mu_c(P[0], phi_p, P[1], P[2])
				return [eq1, eq2, eq3]
			root = fsolve(line_mu, [ps1, ps2, pp2], args=(pp1))
			error = line_mu(root, pp1)
			p1 = np.array([ps1, root[0]])
			p2 = np.array([root[1], root[2]])
			if (np.abs(error) < 1e-6).all(): 
				ax.plot([ps1, root[1]], [1-ps1-root[0],1-root[1]-root[2]], [root[0], root[2]], c="coral", markersize=0.5)
			else:
				bad_line += 1
		elif y_class_pred[i] == 2:
			net_triangle += 1
			ps1 = y_comp_pred[i][0]
			pp1 = y_comp_pred[i][1]
			ps2 = y_comp_pred[i][2]
			pp2 = y_comp_pred[i][3]
			ps3 = y_comp_pred[i][4]
			pp3 = y_comp_pred[i][5]
			root = fsolve(triangle_mu, [ps1, pp1, ps2, pp2, ps3, pp3])
			error = triangle_mu(root)
			# print(f"triangle error = {triangle_mu(root)}", flush=True)
			p1 = np.array([root[0], root[1]])
			p2 = np.array([root[2], root[3]])
			p3 = np.array([root[4], root[5]])
			if (np.abs(error) < 1e-6).all():
				ax.plot([p1[0], p2[0], p3[0], p1[0]], [1-p1[0]-p1[1], 1-p2[0]-p2[1], 1-p3[0]-p3[1], 1-p1[0]-p1[1]], [p1[1], p2[1], p3[1], p1[1]], c="steelblue", markersize=0.5)
			else:
				bad_triangle += 1
		else:
			print(f"y_class_pred[i] = {y_class_pred[i]}", flush=True)
			print(f"Exiting...")
			exit()
	
	if net_line != 0:
		print(f"Fraction of bad line solves = {bad_line/net_line}")
	if net_triangle != 0:
		print(f"Fraction of bad triangle solves = {bad_triangle/net_triangle}")
	fig.savefig(args.o, dpi=1200, bbox_inches="tight")






