import numpy as np
import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# get the struct right
def unfuck_polymer(polymer):
	unfucked_polymer = copy.copy(polymer) # np.asarray([polymer[0,:]])

	for i in range ( polymer.shape[0]-1 ) :
		diff = polymer[i+1,:] - polymer[i,:]

		for j in range(3):
			diff[j] = diff[j] if np.abs(diff[j])==1 else -1*np.sign(diff[j])

		unfucked_polymer[i+1] = unfucked_polymer[i]+diff # np.vstack( (unfucked_polymer, unfucked_polymer[i]+diff ) )

	return unfucked_polymer


# Function to plot a cube
def plot_cube(ax, center, edge_length=1):

	# Create a list of vertices for a cube
	r = edge_length / 2
	vertices = np.array([
		[1, 1, 1],
		[1, 1, -1],
		[1, -1, -1],
		[1, -1, 1],
		[-1, 1, 1],
		[-1, 1, -1],
		[-1, -1, -1],
		[-1, -1, 1]
	]) * r

	vertices += center

	faces = [
		[vertices[0], vertices[1], vertices[2], vertices[3]],
		[vertices[4], vertices[5], vertices[6], vertices[7]],
		[vertices[0], vertices[1], vertices[5], vertices[4]],
		[vertices[2], vertices[3], vertices[7], vertices[6]],
		[vertices[1], vertices[2], vertices[6], vertices[5]],
		[vertices[4], vertices[7], vertices[3], vertices[0]]
	]

	poly3d = Poly3DCollection(faces, facecolors='steelblue', linewidths=1, edgecolors='k', alpha=.25)
	ax.add_collection3d(poly3d)

	return

def filter(A, B):

	dtype = {'names': ['f{}'.format(i) for i in range(A.shape[1])],
			'formats': [A.dtype] * A.shape[1]}
	A_struct = A.view(dtype)
	B_struct = B.view(dtype)

	# Find rows in A that are not in B
	mask = ~np.in1d(A_struct, B_struct)

	# Filter A using the mask
	filtered_A = A[mask]

	return filtered_A

if __name__=="__main__":

	# set size of lattice 
	M = 34

	# set number of blocks
	Nblock = 27

	# these are your potential neighbors
	directions = np.array([ [1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1], \
	[1,1,0], [1,0,1], [1,-1,0], [1,0,-1], [-1,1,0], [-1,0,1], [-1,-1,0], [-1,0,-1], \
	[0,1,1], [0,1,-1], [0,-1,1], [0,-1,-1], \
	[1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1], [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1]], dtype=int)

	# instantiate a random location 
	p0 = [np.random.randint(0, M-1), np.random.randint(0, M-1), np.random.randint(0, M-1)]
	all_particles = np.array([p0], dtype=int)

	# create a store of max_contacts
	max_contact_store = [0]

	# start placing the blocks
	for n in range(Nblock-1):
		all_neighbors = all_particles[:, np.newaxis, :] + directions[np.newaxis, :, :]
		all_neighbors = all_neighbors.reshape(all_particles.shape[0] * directions.shape[0], 3) % M
		all_neighbors = filter(all_neighbors, all_particles)

		max_contacts_rsum = -1
		for idx, neighbor in enumerate(all_neighbors):
			# print((all_particles - neighbor.reshape(1, -1),))
			contacts = np.all(np.isin(all_particles - neighbor.reshape(1, -1) % M, [-1,0,1]), axis=1)
			if len(contacts) > max_contacts_rsum:
				max_contacts_rsum = len(contacts)
				max_idx      = idx

		max_contact_store.append(max_contact_store[-1]+max_contacts_rsum)
		print(f"delta mm = {max_contacts_rsum}")
		all_particles = np.vstack((all_particles, all_neighbors[max_idx]))

	# print(all_particles)
	print(f"Maximum number of contacts = {max_contact_store}", flush=True)

	all_particles = unfuck_polymer(all_particles)

	fig = plt.figure()
	ax  = fig.add_subplot(111, projection='3d')

	# plot each cube
	for particle in all_particles:
		plot_cube(ax, particle)

	# Set the limits of the plot
	ax.set_xlim([np.min(all_particles[:,0])-5, np.max(all_particles[:,0])+5])
	ax.set_ylim([np.min(all_particles[:,1])-5, np.max(all_particles[:,1])+5])
	ax.set_zlim([np.min(all_particles[:,2])-5, np.max(all_particles[:,2])+5])

	# Labels
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	fig.show()
	# fig.savefig("maximal_contacts", dpi=1200, bbox_inches="tight")