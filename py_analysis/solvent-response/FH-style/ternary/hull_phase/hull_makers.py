import numpy as np
import re

if __name__=="__main__":

	P = []
	f = open("db2/files", 'r')
	for line in f:
		match_sc = re.search(r"chisc_([-+]?\d+\.\d+)", line)
		match_ps = re.search(r"chips_([-+]?\d+\.\d+)", line)
		match_pc = re.search(r"chipc_([-+]?\d+\.\d+)", line)
		match_vs = re.search(r"vs_([-+]?\d+\.\d+)", line)
		match_vp = re.search(r"vc_([-+]?\d+\.\d+)", line)
		match_vc = re.search(r"vp_([-+]?\d+\.\d+)", line)
		params = [float(match_sc.group(1)), float(match_ps.group(1)), float(match_pc.group(1)), float(match_vs.group(1)), float(match_vp.group(1)), float(match_vc.group(1))]
		P.append(params)

	P = np.array(P)

	npert = 6
	my_f = open("my_params.txt", 'w')
	for p in P:
		new_params = p + 0.1*np.random.randn(npert, 6)
		for i in range(npert):
			my_f.write(f"vs={new_params[i][0]}\nvc={new_params[i][1]}\nvp={new_params[i][2]}\nchisc={new_params[i][3]}\nchips={new_params[i][4]}\nchipc={new_params[i][5]}\n")
			my_f.write('srun --ntasks=1 --nodes=1 --cpus-per-task=1 --mem-per-cpu=6G --exclusive python ../../compute.py -vs $vs -vc $vc -vp $vp --chisc $chisc --chips $chips --chipc $chipc --mesh 5001 --hull chisc_${chisc}-chips_${chips}-chipc_${chipc}-vs_${vs}-vp_${vp}-vc_${vc}.hull.pkl > chisc_${chisc}-chips_${chips}-chipc_${chipc}-vs_${vs}-vp_${vp}-vc_${vc}.out 2>&1 &\n\n')

