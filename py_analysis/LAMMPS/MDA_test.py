import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis

if __name__=="__main__":
	u = mda.Universe(PSF, DCD)
	analysis = SASAAnalysis(u)
	analysis.run()
	print(analysis.results.total_area)
