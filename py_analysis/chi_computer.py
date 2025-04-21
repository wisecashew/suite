#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np
import re
import argparse

parser = argparse.ArgumentParser (description="Print out chis.")
parser.add_argument ("--path", dest="p", type=str, action='store', help="Path to geometry and energetics file.")
args = parser.parse_args()


def get_info(topology):

	f = open (topology, 'r')
	kT      = "kT"
	Emm_a   = "Emm_a"
	Emm_n   = "Emm_n"
	Ems1_a  = "Ems1_a"
	Ems1_n  = "Ems1_n"
	Ems2_a  = "Ems2_a"
	Ems2_n  = "Ems2_n"
	Es1s2_a = "Es1s2_a"
	Es1s2_n = "Es1s2_n"
	m1m1    = "m1-m1"
	m1s1    = "m1-s1"
	m1s2    = "m1-s2"
	s1s2    = "s1-s2"
	num_re  = "\s+-[0-9]+\.[0-9]+|\s+[0-9]+\.[0-9]+|\s+-[0-9]+\.|\s+[0-9]+\.|\s+-[0-9]+|\s+[0-9]+"
	e_dict = dict()

	for line in f:
		if re.findall (kT, line):
			r = re.findall(num_re, line)
			e_dict["T"] = float(r[0])
		if re.findall (Emm_a, line):
			r = re.findall (num_re, line)
			mm_a = float ( r[0] )
			e_dict["Emm_a"] = float( r[0] )
		elif re.findall ( Emm_n, line):
			r = re.findall (num_re, line)
			e_dict["Emm_n"] = float( r[0] )
		elif re.findall ( Ems1_a, line):
			r = re.findall (num_re, line)
			ms1_a = float( r[0] )
			e_dict["Ems1_a"] = float( r[0] )
		elif re.findall ( Ems1_n, line):
			r = re.findall (num_re, line)
			ms1_n = float( r[0] )
			e_dict["Ems1_n"] = float( r[0] )
		elif re.findall ( Ems2_a, line):
			r = re.findall (num_re, line)
			e_dict["Ems2_a"] = float( r[0] )
		elif re.findall ( Ems2_n, line):
			r = re.findall (num_re, line)
			e_dict["Ems2_n"] = float( r[0] )
		elif re.findall ( Es1s2_a, line):
			r = re.findall (num_re, line)
			e_dict["Es1s2_a"] = float( r[0] )
		elif re.findall ( Es1s2_n, line):
			r = re.findall (num_re, line)
			s1s2_n = float(r[0] )
			e_dict["Es1s2_n"] = s1s2_n
		elif re.findall(m1m1, line):
			typ = re.search(r":(.*)", line)
			e_dict["mm_type"] = typ.group(1)
			typ = typ.group(1)
			if typ == "isotropic":
				e_dict["pv_mm"] = 0
				e_dict["pw_mm"] = 0
			elif typ == "parallel":
				e_dict["pv_mm"] = 1.0
				e_dict["pw_mm"] = 0.25
			elif typ == "antiparallel":
				e_dict["pv_mm"] = 0.35
				e_dict["pw_mm"] = 0.25
		elif re.findall(m1s1, line):
			typ = re.search(r":(.*)", line)
			e_dict["ms1_type"] = typ.group(1)
			typ = typ.group(1)
			if typ == "isotropic":
				e_dict["pv_ms1"] = 0
				e_dict["pw_ms1"] = 0
			elif typ == "parallel":
				e_dict["pv_ms1"] = 1.0
				e_dict["pw_ms1"] = 0.25
			elif typ == "antiparallel":
				e_dict["pv_ms1"] = 0.35
				e_dict["pw_ms1"] = 0.25
			
		elif re.findall(m1s2, line):
			typ = re.search(r":(.*)", line)
			e_dict["ms2_type"] = typ.group(1)
			typ = typ.group(1)
			if typ == "isotropic":
				e_dict["pv_ms2"] = 0
				e_dict["pw_ms2"] = 0
			elif typ == "parallel":
				e_dict["pv_ms2"] = 1.0
				e_dict["pw_ms2"] = 0.25
			elif typ == "antiparallel":
				e_dict["pv_ms2"] = 0.35
				e_dict["pw_ms2"] = 0.25
		elif re.findall(s1s2, line):
			typ = re.search(r":(.*)", line)
			e_dict["s1s2_type"] = typ.group(1)
			typ = typ.group(1)
			if typ == "isotropic":
				e_dict["pv_s1s2"] = 0
				e_dict["pw_s1s2"] = 0
			elif typ == "parallel":
				e_dict["pv_s1s2"] = 1.0
				e_dict["pw_s1s2"] = 0.25
			elif typ == "antiparallel":
				e_dict["pv_s1s2"] = 0.35
				e_dict["pw_s1s2"] = 0.25


	f.close()
	return e_dict

Z     = 26
zmm   = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128) + (1-pw)*np.exp (-1/T * emmn, dtype=np.float128)
zms   = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128) + (1-pw)*np.exp (-1/T * emsn, dtype=np.float128)
zss   = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128) + (1-pw)*np.exp (-1/T * essn, dtype=np.float128)
fmma  = lambda emma, emmn, pw, T: pw*np.exp (-1/T * emma, dtype=np.float128)/zmm(emma, emmn, pw, T)
fmsa  = lambda emsa, emsn, pw, T: pw*np.exp (-1/T * emsa, dtype=np.float128)/zms(emsa, emsn, pw, T)
fssa  = lambda essa, essn, pw, T: pw*np.exp (-1/T * essa, dtype=np.float128)/zss(essa, essn, pw, T)

def compute_chi (e_dict):

	w_mm = e_dict["pv_mm"] * (fmsa(e_dict["Emm_a"], e_dict["Emm_n"], e_dict["pw_mm"], e_dict["T"])*e_dict["Emm_a"] + (1-fmsa(e_dict["Emm_a"], e_dict["Emm_n"], e_dict["pw_mm"], e_dict["T"]))*e_dict["Emm_n"]) + (1-e_dict["pv_mm"]) * (e_dict["Emm_n"])

	w_ms1 = e_dict["pv_ms1"] * (fmsa(e_dict["Ems1_a"], e_dict["Ems1_n"], e_dict["pw_ms1"], e_dict["T"])*e_dict["Ems1_a"] + (1-fmsa(e_dict["Ems1_a"], e_dict["Ems1_n"], e_dict["pw_ms1"], e_dict["T"]))*e_dict["Ems1_n"]) + (1-e_dict["pv_ms1"]) * (e_dict["Ems1_n"])

	w_ms2 = e_dict["pv_ms2"] * (fmsa(e_dict["Ems2_a"], e_dict["Ems2_n"], e_dict["pw_ms2"], e_dict["T"])*e_dict["Ems2_a"] + (1-fmsa(e_dict["Ems2_a"], e_dict["Ems2_n"], e_dict["pw_ms1"], e_dict["T"]))*e_dict["Ems2_n"]) + (1-e_dict["pv_ms2"]) * (e_dict["Ems2_n"])

	w_s1s2 = e_dict["pv_s1s2"] * (fmsa(e_dict["Es1s2_a"], e_dict["Es1s2_n"], e_dict["pw_s1s2"], e_dict["T"])*e_dict["Es1s2_a"] + (1-fmsa(e_dict["Es1s2_a"], e_dict["Es1s2_n"], e_dict["pw_s1s2"], e_dict["T"]))*e_dict["Es1s2_n"]) + (1-e_dict["pv_s1s2"]) * (e_dict["Es1s2_n"])

	# print(f"w_mm  = {w_mm}")
	# print(f"w_ms1 = {w_ms1}")
	# print(f"w_ms2 = {w_ms2}")

	chi_ms1  = (Z-2)/e_dict["T"] * (w_ms1 - 1/2 * w_mm)
	chi_ms2  = (Z-2)/e_dict["T"] * (w_ms2 - 1/2 * w_mm)
	chi_s1s2 = (Z-2)/e_dict["T"] * (w_s1s2)

	return [chi_ms1, chi_ms2, chi_s1s2]

if __name__=="__main__":

	e_dict = get_info(args.p)
	chis   = compute_chi(e_dict)

	print(f"chi_ms1  = {chis[0]}.")
	print(f"chi_ms2  = {chis[1]}.")
	print(f"chi_s1s2 = {chis[2]}.")
