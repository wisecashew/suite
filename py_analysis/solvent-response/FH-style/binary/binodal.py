import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root

if __name__=="__main__":

	mu_p      = lambda phi, m, de, T: np.log(phi) + (1-m)*(1-phi) + de/T*(1-phi)**2*m
	mu_s      = lambda phi, m, de, T: np.log(1-phi) + (1-1/m)*phi + de/T*(phi)**2

	phi     = lambda m,chi: [-(-1 + m-2*chi*m-np.sqrt(-8*chi*m+(1-m+2*chi*m)**2))/(4*chi*m),\
	 -(-1 + m-2*chi*m+np.sqrt(-8*chi*m+(1-m+2*chi*m)**2))/(4*chi*m)]
	chi     = lambda dE, T: dE/T
	T_range = np.logspace(-3,2,10000)

	fig = plt.figure(num=0, figsize=(3,3))
	ax  = plt.axes()
	M   = 32
	dE  = 10 
	spinodal = phi(M, chi(dE, T_range))

	mask0 = np.logical_and(spinodal[0]>0, spinodal[0]<1)
	spinodal[0] = spinodal[0][mask0]
	T0          = T_range[mask0]
	
	mask1 = np.logical_and(spinodal[1]>0, spinodal[1]<1)
	spinodal[1] = spinodal[1][mask1]
	T1          = T_range[mask1]

	phi_c = (-1+np.sqrt(M))/(-1+M)
	T_c   = dE/((np.sqrt(M)+1)**2/(2*M))

	d_dmus_dp1 = lambda x1, m,  de, T:  np.float64(1 - 1/m - 1/(1-x1) + 2*de*x1/T)
	d_dmus_dp2 = lambda x2, m,  de, T:  np.float64(-1 + 1/m + 1/(1-x2) - 2*de*x2/T)
	d_dmus_dT  = lambda x1, x2, m,  de, T: np.float64(-de*(x1)**2/T**2 + de*(x2)**2/T**2)

	d_dmup_dp1 = lambda x1, m,  de, T: np.float64(-1 + m + 1/x1 - 2*m*de*(1-x1)/T)
	d_dmup_dp2 = lambda x2, m,  de, T: np.float64(1 - m - 1/x2 + 2*m*de*(1-x2)/T)
	d_dmup_dT  = lambda x1, x2, m,  de, T: np.float64(-de*m*(1-x1)**2/T**2 + de*m*(1-x2)**2/T**2)

	delta_T  = -0.01 # has to be negative!

	# delta_x1 = 0.01
	# run a cycle until binodal is covered
	ncycles = int(T_c//-delta_T) + 1
	# preloaded steps
	T         = [T_c]
	bin_arm_1 = [phi_c]
	bin_arm_2 = [phi_c]
	# print(ncycles)
	scale = 1
	for i in range(ncycles):
		if T[-1] + delta_T < 0.001:
			break
		print (f"@ i = {i}/{ncycles}...")

		if i > 0:
			# update scheme when outside of crit point
			a3 = -d_dmus_dT (bin_arm_1[-1], bin_arm_2[-1], M, dE, T[-1])*delta_T
			b3 = -d_dmup_dT (bin_arm_1[-1], bin_arm_2[-1], M, dE, T[-1])*delta_T
			a2 = d_dmus_dp2(bin_arm_2[-1], M, dE, T[-1])
			b2 = d_dmup_dp2(bin_arm_2[-1], M, dE, T[-1])
			a1 = d_dmus_dp1(bin_arm_1[-1], M, dE, T[-1])
			b1 = d_dmup_dp1(bin_arm_1[-1], M, dE, T[-1])
			delta_x1 = (a3*b2-a2*b3)/(a1*b2-a2*b1)
			delta_x2 = (a3*b1-a1*b3)/(a2*b1-a1*b2)

		else:
			print(f"i={i}")
			delta_x1 = 0.01
			delta_x2 = -0.01

		print(f"delta_x1 = {delta_x1}")
		print(f"delta_x2 = {delta_x2}")

		def dmu(phi):
			delta_p = (mu_p(phi[0], M, dE, T[-1]+delta_T) - mu_p(phi[1], M, dE, T[-1]+delta_T))/(phi[0]-phi[1])
			delta_s = (mu_s(phi[0], M, dE, T[-1]+delta_T) - mu_s(phi[1], M, dE, T[-1]+delta_T))/(phi[0]-phi[1])
			return [delta_p, delta_s]

		root = fsolve(dmu, [bin_arm_1[-1]+delta_x1, bin_arm_2[-1]+delta_x2], xtol=1e-9)
		if root[0] > 1 or root[0] < 0 or root[1] > 1 or root[1] < 0:
			print("Breaking out...")
			break

		# root = fsolve(dmu, [phi_c+delta_x1, phi_c+delta_x2], xtol=1e-9)
		print(f"p1 = {root[0]} and p2 = {root[1]}")
		print(f"Delta polymer potentials are: {abs(mu_p(root[0], M, dE, T[-1]+delta_T) - mu_p(root[1], M, dE, T[-1]+delta_T))}")
		print(f"Delta solvent potentials are: {abs(mu_s(root[0], M, dE, T[-1]+delta_T) - mu_s(root[1], M, dE, T[-1]+delta_T))}")

		if abs(mu_p(root[0], M, dE, T[-1]+delta_T) - mu_p(root[1], M, dE, T[-1]+delta_T)) > 1e-9:
			continue
		if abs(mu_s(root[0], M, dE, T[-1]+delta_T) - mu_s(root[1], M, dE, T[-1]+delta_T)) > 1e-9:
			continue

		T.append(T[-1]+delta_T)
		bin_arm_1.append(root[0])
		bin_arm_2.append(root[1])



	ax.scatter(spinodal[0], T0, c='coral', s=2)
	ax.scatter(spinodal[1], T1, c='steelblue', s=2)
	ax.scatter(phi_c, T_c, c='darkred', s=2)
	ax.scatter(bin_arm_1, T, c='slategray', s=1)
	ax.scatter(bin_arm_2, T, c='slategray', s=1)
	fig.savefig('spinodal', dpi=1200, bbox_inches='tight')

