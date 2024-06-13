import numpy as np
import mu
import ternary
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
	
	def tangent_tracing_unity(self, BINODALS, uidx):
		
		idx = BINODALS["groupings"][uidx]["alpha"]["idx"]
		print(f"Get the binodal curve...", flush=True, end=' ')

		BINODALS["groupings"][uidx]["alpha"]["binodals"] = [np.empty((0,2)), np.empty((0,2))]
		normal = BINODALS["crit_info"][BINODALS["groupings"][uidx]["alpha"]["idx"]]["norm_vec"]

		# start with mu_ps
		print("Binodal run for ps along s2...")
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx])

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"phi_1 = {phi_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_2 = {phi_2}")
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 
			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!", flush=True)
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
		
		# do the other run 
		print("Binodal run for ps along p2...")
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx])

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))

		# move on to mu_pc
		print("Binodal run for pc along c2...")
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
			
		print("Binodal run for pc along p2...")
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx])
		
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T			

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
			
		print("Binodal run for sc along c2...")
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
		
		print("Binodal run for sc along s2...")
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 
			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))

		along_normal = True
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])

			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 
			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
			

		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx], along_normal)
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T
		
		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))

		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_c1-phi_s1]).T
		phi_2 = np.array([phi_s2, 1-phi_c2-phi_s2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
			
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 
			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
			
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			if (M1 or M2):
				pass 

			elif np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
		
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			M1 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M2 = (ternary.stab_crit(phi_1[:,0], phi_1[:,1], self.vs, self.vc, self.vp, self.chi_ps, self.chi_pc, self.chi_sc) < 0).any() 
			M  = np.logical_or(M1, M2)
			phi_1 = phi_1[~M]
			phi_2 = phi_2[~M]			

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], phi_1))
		
		return
		
	def tangent_tracing_dyad(self, BINODALS, uidx):
		indices = BINODALS["groupings"][uidx]["raw_list"]
		axis    = (self.crits[indices[0]] - self.crits[indices[1]])/np.linalg.norm(self.crits[indices[0]] - self.crits[indices[1]])
		c_crit  = self.crits[indices[1]]
		BINODALS["groupings"][uidx]["alpha"]["binodals"] = [np.empty((0,3)), np.empty((0,3))]
		for count,idx in enumerate(indices):
			
			print(f"Binodal run for ps along s2.", flush=True)
			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx])
			phi_1 = np.array([phi_s1, phi_p1, 1-phi_s1-phi_p1]).T
			phi_2 = np.array([phi_s2, phi_p2, 1-phi_s2-phi_p2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))

			print(f"Binodal run for ps along p2.", flush=True)
			phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx])
			phi_1 = np.array([phi_s1, phi_p1, 1-phi_s1-phi_p1]).T
			phi_2 = np.array([phi_s2, phi_p2, 1-phi_s2-phi_p2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))
			
			print(f"Binodal run for sc along s2.", flush=True)
			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx])
			phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1, phi_c1]).T
			phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2, phi_c2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))
			
			print(f"Binodal run for sc along c2.", flush=True)
			phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx])
			phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1, phi_c1]).T
			phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2, phi_c2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))
			
			print(f"Binodal run for pc along c2.", flush=True)
			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx])
			phi_1 = np.array([1-phi_p1-phi_c1, phi_p1, phi_c1]).T
			phi_2 = np.array([1-phi_p2-phi_c2, phi_p2, phi_c2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))

			print(f"Binodal run for pc along p2.", flush=True)
			phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx])
			phi_1 = np.array([1-phi_p1-phi_c1, phi_p1, phi_c1]).T
			phi_2 = np.array([1-phi_p2-phi_c2, phi_p2, phi_c2]).T

			adj_phi_1 = (phi_1[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_1[:,0:2]-c_crit[0:2]).reshape(-1,1)
			adj_phi_2 = (phi_2[:,0:2]-c_crit[0:2])/np.linalg.norm(phi_2[:,0:2]-c_crit[0:2]).reshape(-1,1)

			mask1 = (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == -1)
			mask2 = (np.sign(np.cross(axis[0:2], adj_phi_2[:,0:2])) == 1) & (np.sign(np.cross(axis[0:2], adj_phi_1[:,0:2])) == -1)
			
			pos_phi = np.vstack((phi_1[mask1], phi_2[mask2]))
			neg_phi = np.vstack((phi_2[mask1], phi_1[mask2]))

			BINODALS["groupings"][uidx]["alpha"]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][0], neg_phi))
			BINODALS["groupings"][uidx]["alpha"]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx]["alpha"]["binodals"][1], pos_phi))
		
		return

	def tangent_tracing_triumvirate(self, BINODALS, probe, uidx):
		
		idx = BINODALS["groupings"][uidx][probe]["idx"]
		print(f"Get the binodal curve...", flush=True, end=' ')

		BINODALS["groupings"][uidx][probe]["binodals"] = [np.empty((0,2)), np.empty((0,2))]
		normal = BINODALS["crit_info"][BINODALS["groupings"][uidx][probe]["idx"]]["norm_vec"]
		print(f"crit   = {self.crits[idx]}")
		print(f"normal = {normal}")

		# start with mu_ps
		print("Binodal run for ps along s2...")
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx])

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"phi_1 = {phi_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!", flush=True)
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
		
		# do the other run 
		print("Binodal run for ps along p2...")
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx])

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))

		# move on to mu_pc
		print("Binodal run for pc along c2...")
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
			
		print("Binodal run for pc along p2...")
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx])
		
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T			

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
			
		print("Binodal run for sc along c2...")
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"phi_1 = {phi_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
		
		print("Binodal run for sc along s2...")
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx])
		
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))

		along_normal = True
		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_s2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
			

		phi_s1, phi_s2, phi_p1, phi_p2 = self.sym_mu_ps.binodal_run_in_p2(self.crits[idx], along_normal)
		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, phi_p1]).T
		phi_2 = np.array([phi_s2, phi_p2]).T
		
		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))

		print(f"Get the binodal curve...", flush=True, end=' ')
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_s2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_c1-phi_s1]).T
		phi_2 = np.array([phi_s2, 1-phi_c2-phi_s2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
			
		phi_s1, phi_s2, phi_c1, phi_c2 = self.sym_mu_sc.binodal_run_in_c2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([phi_s1, 1-phi_s1-phi_c1]).T
		phi_2 = np.array([phi_s2, 1-phi_s2-phi_c2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print(f"In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
			
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_p2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
		
		phi_p1, phi_p2, phi_c1, phi_c2 = self.sym_mu_pc.binodal_run_in_c2(self.crits[idx], along_normal)

		print("Found a complete binodal. Moving on...")
		phi_1 = np.array([1-phi_p1-phi_c1, phi_p1]).T
		phi_2 = np.array([1-phi_p2-phi_c2, phi_p2]).T

		if len(phi_1) > 1:
			delta_1 = (self.crits[idx][0:2]-phi_1[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_1[1][0:2])
			delta_2 = (self.crits[idx][0:2]-phi_2[1][0:2])/np.linalg.norm(self.crits[idx][0:2]-phi_2[1][0:2])
			print(f"delta_1 = {delta_1}")
			print(f"delta_2 = {delta_2}")
			print(f"phi_1 = {phi_1}")
			print(f"phi_2 = {phi_2}")

			if np.sign(np.cross(normal, delta_1)) * np.sign(np.cross(normal, delta_2)) == -1:
				print("In!")
				if np.sign(np.cross(normal, delta_1)) == 1:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_1))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_2))
				else:
					BINODALS["groupings"][uidx][probe]["binodals"][0] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][0], phi_2))
					BINODALS["groupings"][uidx][probe]["binodals"][1] = np.vstack((BINODALS["groupings"][uidx][probe]["binodals"][1], phi_1))
		
		return
		