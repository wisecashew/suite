#!/home/satyend/.conda/envs/phase/bin/python

import MDAnalysis as mda
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis

if __name__=="__main__":

	topo = "pnipam_tip4p-ice-water.r1.data"
	traj = "coords.lammpstrj"
	dt = 0.001

	# generate the universe object
	u = mda.Universe(topo, traj, format="LAMMPSDUMP", lengthunit="angstrom", timeunit="fs", dt=dt)
	print(f"Dimensions of the universe: {u.dimensions}", flush=True)

	# time to play with molecules in the universe
	print("Begin selecting some atoms for the polymer...", end=' ', flush=True)
	polymer = u.select_atoms('resid 1')
	print("done!")

	analysis = SASAAnalysis(u)
	analysis.run()

	# set up a store for the radius of gyrationa and SASA
	Rg        = []

	for idx,ts in enumerate(u.trajectory):
		polymer.unwrap()
		Rg.append(polymer.radius_of_gyration(wrap=False))
		if idx % 5 == 0:
			print(f"@timestep {idx}, Rg = {Rg[-1]}.", flush=True)
		

	'''
	# Access and print atoms
	print(f"Atoms:")
	for typ in u.atoms.types:
		print(f"Atom type = {typ}")

	print(f"Number of molecules in the system is: {len(u.residues)}")
	
	# Access and print bonds
	print("Bonds:")
	for bond in u.bonds:
		print(f"Bond between atom {bond.atoms[0].id} and atom {bond.atoms[1].id}")

	# Access and print angles
	print("\nAngles:")
	for angle in u.angles:
		print(f"Angle between atom {angle.atoms[0].id}, atom {angle.atoms[1].id}, and atom {angle.atoms[2].id}")

	# Access and print dihedrals
	print("\nDihedrals:")
	for dihedral in u.dihedrals:
		print(f"Dihedral between atom {dihedral.atoms[0].id}, atom {dihedral.atoms[1].id}, "
			f"atom {dihedral.atoms[2].id}, and atom {dihedral.atoms[3].id}")

	# Access and print impropers
	print("\nImpropers:")
	for improper in u.impropers:
		print(f"Improper between atom {improper.atoms[0].id}, atom {improper.atoms[1].id}, "
			f"atom {improper.atoms[2].id}, and atom {improper.atoms[3].id}")
	'''

