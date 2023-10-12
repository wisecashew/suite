import numpy as np 

def unwrap_polymer(polymer, bonds, box_dims):

	untested = list(range(len(polymer)))
	tested   = []
	queue    = []

	while untested:
		
		wait = []
		if not queue:
			queue.append(untested[0])

		for i in queue:
			neighbors = bonds[i]
			neighbors = [ni for ni in neighbors if ni not in tested]
			ri = polymer[i]
			for j in neighbors:
				rj = polymer[j]
				dr = rj - ri
				shift = np.round(dr/box_dims)
				polymer[j] -= shift*box_dims
			tested.append(i)
			untested.remove(i)
			wait.extend(neighbors)
		queue = list(set(wait[:]))

	return polymer 

if __name__=="__main__":

	polymer           = np.array([[0,0,0],[1,0,0],[9,0,0]], dtype=np.float64) # np.array([[9,1,1],[9,0,1],[10,0,1],[1,0.5,1],[8,1,1],[9.8,2,1],[0.8,2,1],[9.8,3,1],[0.3,4,1],[1.3,4,1]])
	unwrapped_polymer = np.array([[0,0,0],[1,0,0],[-1,0,0]], dtype=np.float64)# np.array([[9,1,1],[9,0,1],[10,0,1],[11,0.5,1],[8,1,1],[9.8,2,1],[10.8,2,1],[9.8,3,1],[10.3,4,1],[11.3,4,1]])

	bonds ={0:[1,2], 1:[0], 2:[0]}
	box_dims = np.array([10,10,10])

	u_polymer = unwrap_polymer(polymer, bonds, box_dims)
	print(u_polymer)


