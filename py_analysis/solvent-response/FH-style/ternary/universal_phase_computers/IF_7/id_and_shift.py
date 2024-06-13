import pickle 
import shutil
import os

if __name__=="__main__":
	
	l_files = os.listdir("stable_islands/.")
	for f in l_files:
		si_filename      = f
		ui_filename      = f.replace(".stable_islands.", ".unstable_islands.")
		si               = open("stable_islands/"+si_filename, 'rb')
		stable_islands   = pickle.load(si)
		ui               = open("unstable_islands/"+ui_filename, 'rb')
		unstable_islands = pickle.load(ui)
		if len(stable_islands) == 1 and len(unstable_islands) == 1:
			shutil.copy("stable_islands/"+f  , "REGULAR2/stable_islands/.")
			f_ = f.replace(".stable_islands.", ".unstable_islands.")
			g_ = f.replace(".stable_islands.", ".crits.")
			shutil.copy("unstable_islands/"+f_, "REGULAR2/unstable_islands/.")
			try:
				shutil.copy("crits/"+g_, "REGULAR2/crits/.")
			except FileNotFoundError:
				print("Crit file unavailable.")
		# elif len(stable_islands) > 1:
		# 	shutil.copy("stable_islands/"+f, "BIGMEM/stable_islands/.")
		# 	f_ = f.replace(".stable_islands.", ".unstable_islands.")
		# 	g_ = f.replace(".stable_islands.", ".crits.")
		# 	shutil.copy("unstable_islands/"+f_, "BIGMEM/unstable_islands/.")
		# 	shutil.copy("crits/"+g_, "BIGMEM/crits/.")
		# else:
		# 	print(f"Somethings wrong with {f}")
		# 	exit()
		si.close()
		ui.close()
