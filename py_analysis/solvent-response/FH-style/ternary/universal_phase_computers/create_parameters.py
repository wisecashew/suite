#!/home/satyend/.conda/envs/phase/bin/python

import numpy as np

if __name__=="__main__":

	chisc = np.array([0, 1, 2, 3, 4, 5])
	chips = np.array([0, 1, 2, 3, 4, 5])
	chipc = np.array([0, 1, 2, 3, 4, 5])

	vs = np.array([1, 2, 3, 4, 5])
	vc = np.array([1, 2, 3, 4, 5])
	vp = np.array([1, 2, 3, 4, 5])

	params = []

	for c_sc in chisc:
		for c_ps in chips:
			for c_pc in chipc:
				for s_ in vs:
					for c_ in vc:
						for p_ in vp:
							select = [c_sc, c_ps, c_pc, s_, c_, p_]
							params.append(select)
	
	idxes  = np.random.choice(range(len(params)), size=1000, replace=False)
	params = np.array(params)
	params = params[idxes]
	p_file = open("params.txt", 'w')
	'''
	p_file.write("for chisc in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][0]}; do")
		else:
			p_file.write(f"{params[i][0]} ")
	p_file.write("\n")

	p_file.write("for chips in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][1]}; do")
		else:
			p_file.write(f"{params[i][1]} ")
	p_file.write("\n")

	p_file.write("for chipc in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][2]}; do")
		else:
			p_file.write(f"{params[i][2]} ")
	p_file.write("\n")

	p_file.write("for vs in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][3]}; do")
		else:
			p_file.write(f"{params[i][3]} ")
	p_file.write("\n")

	p_file.write("for vc in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][4]}; do")
		else:
			p_file.write(f"{params[i][4]} ")
	p_file.write("\n")

	p_file.write("for vp in ")
	for i in range(len(params)):
		if i == len(params)-1:
			p_file.write(f"{params[i][5]}; do")
		else:
			p_file.write(f"{params[i][5]} ")
	p_file.write("\n")
	'''
	for p in params:
		p_file.write(f"srun --ntasks=1 --cpus-per-task=1 --nodes=1 --mem-per-cpu=2GB python cisland.py \\ \n")
		p_file.write(f"--ternary --chisc {p[0]} --chips {p[1]} --chipc {p[2]} -vs {p[3]} -vc {p[4]} -vp {p[5]} \\ \n")
		p_file.write(f"--draw-edges-spinodal --mesh-pkl mesh.pkl \\ \n")
		p_file.write(f"--island-stable-pkl IF_7/stable_islands/chisc_{p[0]}_chips_{p[1]}_chipc_{p[2]}-vs_{p[3]}-vc_{p[4]}-vp_{p[5]}.stable_islands.pkl \\ \n")
		p_file.write(f"--island-unstable-pkl IF_7/unstable_islands/chisc_{p[0]}_chips_{p[1]}_chipc_{p[2]}-vs_{p[3]}-vc_{p[4]}-vp_{p[5]}.unstable_islands.pkl \\ \n")
		p_file.write(f"--search-density 1000 --plot-crits > IF_7/output_chisc_{p[0]}_chips_{p[1]}_chipc_{p[2]}_vs_{p[3]}_vc_{p[4]}_vp_{p[5]}.out 2>&1 & \n\n")

	p_file.write("wait\n")
	p_file.close()
