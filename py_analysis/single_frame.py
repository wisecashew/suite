#!/usr/bin/env python3
import re 
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
import argparse 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


parser = argparse.ArgumentParser(description='Reads a trajectory and generates a movie of all the selected frames in the trajectory.')
parser.add_argument('-x', dest='x', action='store', type=int, help='Length of cell along the x axis.')
parser.add_argument('-y', dest='y', action='store', type=int, help='Length of cell along the y axis.')
parser.add_argument('-z', dest='z', action='store', type=int, help='Length of cell along the z axis.')
parser.add_argument('-p', dest='p', action='store', type=int, help='Extract coordinates given at this move.')
parser.add_argument('-n', metavar='stepX.png', dest='n', action='store', type=str, help='Name of image file to be created.')
parser.add_argument('--traj', metavar='coords.txt', dest='traj', action='store', type=str, help='Name of trajectory file.')
parser.add_argument('-w', dest='w', action='store', type=int, default=0, help='Enter 0 if you do not want the whole box view. Enter 1 if you want the whole box view.')

args = parser.parse_args() 

directions = [ [1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1], \
    [1,1,0], [1,0,1], [1,-1,0], [1,0,-1], [-1,1,0], [-1,0,1], [-1,-1,0], [-1,0,-1], \
    [0,1,1], [0,1,-1], [0,-1,1], [0,-1,-1], \
    [1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1], [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1]]



if (args.w != 0 and args.w !=1 ):
    print("Bad args.w value. Only 0 or 1 allowed.")
    exit()

################################
### functions 
################################

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


def extract_loc_from_string (a_string):
    loc = [int(word) for word in a_string.split() if word.isdigit() ] 
    return np.asarray(loc) 

def unfuck_polymer (polymer, x, y, z): 
    unfucked_polymer = np.asarray([polymer[0,:]])
    
    for i in range ( polymer.shape[0]-1 ) : 
        diff = polymer[i+1,:] - polymer[i,:]
    
        for j in range(3):
            diff[j] = modified_modulo(diff[j], x)
    
        unfucked_polymer = np.vstack( (unfucked_polymer, unfucked_polymer[i]+diff ) ) 
    
    return unfucked_polymer  

def modified_modulo(divident, divisor):
    midway = divisor/2
    if (divident%divisor > midway):
        result = (divident%divisor)-divisor
        return result
    else:
        return divident%divisor


def get_ne_list (unfucked_polymer, x, y, z):

	ne_list = [] 
	for loc in unfucked_polymer: 
		for drtn in directions: 
			ne = np.asarray (loc) + np.asarray (drtn) 
			ne[0] = ne[0]%x 
			ne[1] = ne[0]%y 
			ne[2] = ne[2]%z 
			ne_list.append(ne) 

	# make sure list is unique 
	unique_ne_list = [ ] 

	for ne_ in ne_list:
		if len(unique_ne_list) == 0:
			unique_ne_list.append(ne)
		else:
			_var = (np.linalg.norm((np.array(unique_ne_list) - ne_), axis=1) < 1e-4).any()
			if _var:
				continue
			else:
				unique_ne_list.append(ne) 

	return unique_ne_list 



def get_solvation_shell (ne_list, unfucked_polymer): 
    
    for loc in unfucked_polymer:
        if ne_list.count(loc) > 0:
            ne_list.remove(loc) 

    return ne_list 

################################
################################



st_b_str = "Dumping coordinates at step"
pmer_num_str = "Dumping coordinates of Polymer #" 
start_str = "START" 
end_str_1 = "END"
end_str_2 = "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n"



# master_dict will be the dictionary which will contain polymer coordinates at each step 
# this will be a large data structure 

master_dict = {} 

step_flag = 0 
pmer_flag = -1 
end_step_flag = 0
step_num = 0 


coord_file = open( args.traj, 'r' )
coord_file = coord_file.readlines()  

for line in coord_file:
    if ( re.search(st_b_str, line)):
        
        step_num = int( (extract_loc_from_string (line.replace('.', ' ') ) ) )
        # print(step_num)
        master_dict[step_num] = {} 
        
        step_flag = 1
        pmer_flag = -1
        end_step_flag = 0
        continue 
    
    elif  (re.search (start_str, line)  ):
        continue 
    
    elif ( re.search (pmer_num_str, line )):
        # print("hit.")
        pmer_flag += 1
        master_dict[step_num][pmer_flag] = np.empty ( (0,3) )  
        continue 
    
    elif ( re.search(end_str_1, line ) ):
        continue
        
    elif ( re.search( end_str_2, line) ):
        end_step_flag = 1
        step_flag = 0
        pmer_flag = -1 
        continue
        
    else:
        # print(pmer_flag)
        monomer_coords = extract_loc_from_string ( line ) 
        master_dict[step_num][pmer_flag] = np.vstack( (master_dict[step_num][pmer_flag], monomer_coords[0:-1]) )
        continue 



fig=plt.figure(figsize=(4,4)) 
ax = fig.add_subplot(111, projection='3d')

# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.set_zticklabels([]) 


xmin, ymin, zmin = 10000, 10000, 10000
xmax, ymax, zmax = -10000, -10000, -10000

step_to_extract = args.p

for key in master_dict[step_to_extract]:
    deg_poly = np.shape( master_dict[step_to_extract][key] )[0]
    x_coords = []
    y_coords = []
    z_coords = []
    unfucked_polymer = unfuck_polymer( master_dict[step_to_extract][key], args.x, args.y, args.z )
    for i in range(deg_poly):
        
        x_coords.append( unfucked_polymer[i][0] )
        y_coords.append( unfucked_polymer[i][1] )
        z_coords.append( unfucked_polymer[i][2] )
    
    if (np.min(x_coords) < xmin):
        xmin = np.min(x_coords) 
    if (np.min(y_coords) < ymin):
        ymin = np.min(y_coords) 
    if (np.min(z_coords) < zmin):
        zmin = np.min(z_coords) 
    
    if (np.max(x_coords) > xmax):
        xmax = np.max(x_coords)
    if (np.max(y_coords) > ymax):
        ymax = np.max(y_coords) 
    if (np.max(z_coords) > zmax):
        zmax = np.max(z_coords)
        
    ax.plot(x_coords, y_coords, z_coords, c='C1')
    for i in range(len(x_coords)):
        center = np.array([x_coords[i], y_coords[i], z_coords[i]])
        plot_cube(ax, center)

# get the solvation shell... 
# ne_list = get_ne_list ( unfucked_polymer, args.x, args.y, args.z )
# solvation_shell = get_solvation_shell ( ne_list, unfucked_polymer )
# ax.scatter ( solvation_shell[:][0], solvation_shell[:][1], solvation_shell[:][2], marker='o', c='g', alpha=0.5 )

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
if (args.w==0):
    ax.set_xlim3d(left= xmin-1, right=xmax+1)
    ax.set_ylim3d(bottom= ymin-1, top=ymax+1)
    ax.set_zbound(lower= zmin-1, upper =zmax+1)
else:
    ax.set_xlim3d(left= -args.x, right=args.x)
    ax.set_ylim3d(bottom= -args.y, top=args.y)
    ax.set_zbound(lower= -args.z, upper=args.z)

plt.grid(False)
plt.savefig(args.n, dpi=1200)
plt.show()
    
