#!/usr/licensed/anaconda3/2020.7/bin/python

# PPG stands for Polymer Positions Generator. 
# given the number of polymers and dimensions of the box, PPG should be able to create valid 
# simulation boxes 
 
import numpy as np 
import argparse 
import copy 
import sys 
import re

sys.setrecursionlimit(100000)

# use argparse to set up the options 

# set up the description 
parser = argparse.ArgumentParser(description="Get cubical box dimensions (x, y, z) and number of polymers with degree of polymerization to get a valid location file. ")

# -x is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-x", dest='x', metavar="X dimension", type=int, help="Enter length of x-edge of box") 

# -y is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-y", dest='y', metavar="Y dimension", type=int, help="Enter length of y-edge of box") 

# -z is going to be your flag for inputting the .pdb file you want to convert 
parser.add_argument("-z", dest='z', metavar="Z dimension", type=int, help="Enter length of z-edge of box") 

# -d is the degree of polymerization of the polymers in the box 
parser.add_argument("-d", dest='d', metavar="Degree of polymerization", type=int, help="Enter the degree of polymerization of the monodisperse polymers in box")

# -n is the number of polymers in the box 
parser.add_argument("-n", dest='n', metavar="Number of polymers in the box", type=int, help="Enter the number of polymers in the box")

# -o is going to be your flag for defining your output file 
parser.add_argument("-o", dest='o', metavar="positions.txt", type=str, help="Enter name of file which will be filled with coordinates") 

parser.add_argument("--arch", dest='arch', metavar="A#B#", type=str, help="Enter type of architecture for copolymer." ) 

args = parser.parse_args() 

# print("The file path you have provided is " + args.i[0]) 
# print("The name of your output file is " + args.o[0]) 

directions = np.asarray([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]] ) 


##########################################################
### create a lattice 
##########################################################

def OccupancyMap(x, y, z): 
    OccupancyMap = { }    
    for i in range(x):
        for j in range(y):
            for k in range(z): 
                lpoint = np.asarray( [int(i), int(j), int(k)] )
                OccupancyMap [lpoint] = 0 
    
    return OccupancyMap


##########################################################    

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

def sarw_add(existing, pmer, DOP, x, y, z):
    
    if (DOP == 0): 
        return existing, pmer 
    
    np.random.shuffle(directions) 
    # print("existing is:")
    # print(existing)
    
    for increment in directions: 
        possible_next = existing[-1] + increment 
        possible_next = impose_pbc (possible_next, x, y, z) 
        
        if (not (possible_next.tolist() in existing.tolist() ) ):
            DOP -= 1 
            # print(pmer)
            pmer = np.vstack ( (pmer, possible_next) )
            existing = np.vstack( (existing, possible_next) )
            
            return sarw_add (existing, pmer, DOP, x, y, z) 
    
    print("No good solutions found...")
    return existing, pmer
    
##########################################################

# necessary inputs 
if __name__=="__main__":

    deg_poly = args.d
    num_poly = args.n
    xlen = args.x 
    ylen = args.y
    zlen = args.z
    arch = args.arch 

    # seed monomer 
    polymer = np.asarray( [[0,0,0]] ) 

    # define a polymer 
    polymer = sarw_pbc (polymer, deg_poly-1, xlen, ylen, zlen) 

    plocs = copy.deepcopy(polymer) 

    # appending to the global polymer list 
    PolymerList = [polymer]

    # add polymers in box 
    for i in range(num_poly-1):
        plocs, new_polymer = sarw_add (plocs, np.empty((0,3)), deg_poly, xlen, ylen, zlen) 
        PolymerList.append (new_polymer) 

    
    A_info = re.findall("A\d+", arch)
    
    A_n = int(A_info[0][1:])
    B_info = re.findall("B\d+", arch)
    B_n = int(B_info[0][1:])

    # l2m = {'A': '1', 'B': '2'}
    arch_length = A_n+B_n 

    def i2m(idx, al_, A_n_, B_n_):
        if idx%al_ < A_n_:
            return '1'
        else:
            return '2' 


    mtypes = [] 
    for i in range (deg_poly):
        mtypes.append("m"+i2m(i, arch_length, A_n, B_n)) 

    # send locations of polymer to coordinate file 
    
    f = open(args.o, "w")

    for i in range(len(PolymerList)):

        f.write("START POLYMER {:1.0f}\n".format(int(i)) )
        for coords,mtype in zip(PolymerList[i], mtypes):
            f.write(mtype + ": " + str(int(coords[0])) + " " + str(int(coords[1])) + " " + str(int(coords[2])) + "\n")
            
        if (i==len(PolymerList)-1):
            f.write("END POLYMER {:1.0f}".format(int(i)) )
        else:
            f.write("END POLYMER {:1.0f}\n".format (int(i)) )

    f.close()                         
    

