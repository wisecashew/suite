#!/usr/bin/env python3
import re 
import numpy as np
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description='Reads a trajectory and generates a movie of all the selected frames in the trajectory.')
parser.add_argument('-x', dest='x', action='store', type=int, help='Length of cell along the x axis.')
parser.add_argument('-y', dest='y', action='store', type=int, help='Length of cell along the y axis.')
parser.add_argument('-z', dest='z', action='store', type=int, help='Length of cell along the z axis.')
parser.add_argument('-n', dest='n', action='store', type=str, help='Name of mp4 file to be created.')
parser.add_argument('--traj', dest='traj', action='store', type=str, help='Name of trajectory file.')
parser.add_argument('-f1', dest='f1', action='store', type=int, help='Index of first frame.')
parser.add_argument('-f2', dest='f2', action='store', type=int, help='Index of final frame.')
parser.add_argument('-s', dest='s', action='store', type=int, help='Steps between two frames.')

args = parser.parse_args() 

################################
### functions 
################################

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


Writer = animation.writers["ffmpeg"]
writer = Writer(fps=2, metadata=dict(artist='Me'))

fig=plt.figure(figsize=(4,4)) 
ax = fig.add_subplot(111, projection='3d')

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([]) 





def animate(step_to_extract):
    xmin, ymin, zmin = 10000, 10000, 10000
    xmax, ymax, zmax = -10000, -10000, -10000
    ax.clear() 
    for key in master_dict[step_to_extract]:
        # print(key)
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
        
        # print("step number is " + str(step_to_extract))
        # print(key)
        # print(x_coords)    
        ax.plot(x_coords, y_coords, z_coords, c='C1')
        ax.scatter(x_coords, y_coords, z_coords, marker='o', c='g', alpha=0.5)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_xlim3d(left=xmin-1, right=xmax+1)
    ax.set_ylim3d(bottom=ymin-1, top=ymax+1)
    ax.set_zbound(lower=zmin-1, upper =zmax+1)
    ax.set_title("Step number = " + str(step_to_extract))
    im = ax 

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=range(args.f1, args.f2+1,args.s), interval=1000, repeat=False, blit=False)

ani.save(args.n, writer=writer)
plt.show()
    
