# PPG stands for Polymer Positions Generator. 
# given the number of polymers and dimensions of the box, PPG should be able to create valid 
# simulation boxes 

import numpy as np 
import argparse 

# use argparse to set up the options 
'''
# set up the description 
parser = argparse.ArgumentParser(description="Get cubical box dimensions (x, y, z) and number of polymers with degree of polymerization to get a valid location file. ")

# -x is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-x", metavar="X dimension", type=int, nargs=1, help="Enter length of x-edge of box") 

# -y is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-y", metavar="Y dimension", type=int, nargs=1, help="Enter length of y-edge of box") 

# -z is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-z", metavar="Z dimension", type=int, nargs=1, help="Enter length of z-edge of box") 


# -o is going to be your flag for defining your output file 
parser.add_argument("-o", metavar="positions.txt", type=str, nargs=1, help="Enter name of file which will be filled with coordinates") 

# parse the aruments 
# args.i[0] is the string you get to use in this code for the address of your input 
# args.o[0] is the string you get to use in this code for your output 

args = parser.parse_args() 

print("The file path you have provided is " + args.i[0]) 
print("The name of your output file is " + args.o[0]) 
'''
directions = np.asarray([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]] ) 

##########################################################
### impose periodic boundary conditions on monomer
##########################################################
 
def impose_pbc(loc, x, y, z):
    for i in range(3):
        if (i==0):
            loc[i] = (( loc[i]%x ) +x ) %x
        
        elif (i==1):
            loc[i] = (( loc[i]%y ) +y ) %y
            
        else:
            loc[i] = (( loc[i]%z ) +z ) %z
            
    return loc 

##########################################################    


##########################################################
### create self-avoiding polymer 
##########################################################
def sarw_pbc(positions, DOP, x, y, z):
    
    if (DOP==0):
        return positions
    
    np.random.shuffle(directions) 
    for increment in directions:
        possible_next = positions[-1] + increment
        possible_next = impose_pbc(possible_next, x, y, z) 
        # print(possible_next) 
        
        if (not (possible_next.tolist() in positions.tolist() ) ):
            DOP -= 1
            positions = np.vstack( (positions, possible_next) )
            
            return sarw_pbc (positions, DOP, x, y, z)
        
    print("No good solution found...")
    return positions
##########################################################


dop_list = [8, 16, 24, 32, 48, 64]

for dop in dop_list:
    # box dimensions 
    x = dop+2
    y = dop+2 
    z = dop+2

    polymer = np.asarray([[0,0,0]]) # seed location 

    polymer = sarw_pbc(polymer, dop-1, x, y, z)
 
    # print(polymer) 

    # send locations of polymer to coordinate file 
    f = open("positions_"+ str(dop) + ".txt", "w")
    f.write("START POLYMER 1\n")
    for coords in polymer:
        f.write(str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")
    f.write("END POLYMER 1")
    f.close()             
            


    

