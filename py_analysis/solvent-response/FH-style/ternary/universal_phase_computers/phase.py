import numpy as np
import mu
import spinodal

class Phase:
	def __init__(self, inputs, crits=None):
		self.chi_sc    = inputs["chi_sc"]
		self.chi_ps    = inputs["chi_ps"]
		self.chi_pc    = inputs["chi_pc"]
		self.vs        = inputs["vs"]
		self.vc        = inputs["vc"]
		self.vp        = inputs["vp"]
		self.spinodal  = spinodal.Spinodal(inputs)
		self.crits     = self.spinodal.crits
		self.sym_mu_ps = mu.sym_mu_ps(inputs, self.spinodal)
		self.sym_mu_sc = mu.sym_mu_sc(inputs, self.spinodal)
		self.sym_mu_pc = mu.sym_mu_pc(inputs, self.spinodal)
		return
	
	def tangent_tracing(self, BINODALS):
		bin_count = 0
		for idx, c in enumerate(self.crits):
			
			for test_idx in BINODALS["groupings"]:
				if idx in BINODALS["groupings"][test_idx]["raw_list"]:
					uidx = test_idx
					break
			

			print(f"Get the binodal curve...", flush=True, end=' ')

			# start with mu_ps
			print("Binodal run for ps along s2...")
			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx])

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='gold', s=3)
				# ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='pink', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, phi_p1]).T
				phi_2 = np.array([phi_s2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			# do the other run 
			print("Binodal run for ps along p2...")
			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx])

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='gold', s=3)
				# ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='pink', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, phi_p1]).T
				phi_2 = np.array([phi_s2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			
			# move on to mu_pc
			print("Binodal run for pc along c2...")
			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx])
			
			if (phi_p1[-1]<1e-2 or phi_c1[-1]<1e-2 or (1-phi_c1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting 
				# ax.scatter(1-phi_p1-phi_c1, phi_c1, phi_p1, c='gold', s=3)
				# ax.scatter(1-phi_p2-phi_c2, phi_c2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
				phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue



			print("Binodal run for pc along p2...")
			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx])
			
			if (phi_p1[-1]<1e-2 or phi_c1[-1]<1e-2 or (1-phi_c1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting 
				# ax.scatter(1-phi_p1-phi_c1, phi_c1, phi_p1, c='gold', s=3)
				# ax.scatter(1-phi_p2-phi_c2, phi_c2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			
			print("Binodal run for sc along c2...")
			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx])
			
			if (phi_s1[-1]<1e-2 or phi_c1[-1]<1e-2 or (1-phi_c1[-1]-phi_s1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting 
				# ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='gold', s=3)
				# ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
				phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			
			print("Binodal run for sc along s2...")
			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx])
			
			if (phi_s1[-1]<1e-2 or phi_c1[-1]<1e-2 or (1-phi_c1[-1]-phi_s1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting 
				# ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='gold', s=3)
				# ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
				phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

		
			along_normal = True
			print(f"Get the binodal curve...", flush=True, end=' ')
			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='gold', s=3)
				# ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, phi_p1]).T
				phi_2 = np.array([phi_s2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, 1-phi_s1-phi_p1, phi_p1, c='gold', s=3)
				# ax.scatter(phi_s2, 1-phi_s2-phi_p2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, phi_p1]).T
				phi_2 = np.array([phi_s2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue
			
			print(f"Get the binodal curve...", flush=True, end=' ')
			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='gold', s=3)
				# ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, 1-phi_p1-phi_s1]).T
				phi_2 = np.array([phi_s2, 1-phi_p2-phi_s2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2) or (min_dist < 0.01)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(phi_s1, phi_c1, 1-phi_s1-phi_c1, c='gold', s=3)
				# ax.scatter(phi_s2, phi_c2, 1-phi_s2-phi_c2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
				phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue
		
			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(1-phi_p1-phi_c1, phi_c1, phi_p1, c='gold', s=3)
				# ax.scatter(1-phi_p2-phi_c2, phi_c2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
				phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue

			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx], along_normal)

			if (phi_s1[-1]<1e-2 or phi_p1[-1]<1e-2 or (1-phi_s1[-1]-phi_p1[-1]<1e-2) or (min_dist < 0.01)) and len(phi_s1) > 10:
				# start plotting out the curve
				# ax.scatter(1-phi_c1-phi_p1, phi_c1, phi_p1, c='gold', s=3)
				# ax.scatter(1-phi_c2-phi_p2, phi_c2, phi_p2, c='lime', s=3)
				print("Found a complete binodal. Moving on...")
				phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
				phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T
				BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
				continue
			
			BINODALS["groupings"][uidx]["center"]["binodals"] = [phi_1, phi_2]
		
		return 