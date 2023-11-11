import numpy as np
import mu
import spinodal

class Phase:
	def __init__(self, inputs, crits=None):
		self.chi_sc   = inputs["chi_sc"]
		self.chi_ps   = inputs["chi_ps"]
		self.chi_pc   = inputs["chi_pc"]
		self.vs       = inputs["vs"]
		self.vc       = inputs["vc"]
		self.vp       = inputs["vp"]
		self.spinodal = spinodal.Spinodal(inputs)
		self.crits    = self.spinodal.crits
		self.mu_ps    = mu.Mu_PS(inputs)
		self.mu_pc    = mu.Mu_PC(inputs)
		self.mu_sc    = mu.Mu_SC(inputs)
		return

