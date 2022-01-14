import os 
import re 
import matplotlib.pyplot as plt
import numpy as np 

cwd = os.getcwd() 
files = os.listdir() 

T = []
Rg = [] 

for file in files:
    if (re.match('stat_kt', file)):
        T.append( float(re.findall("\d+\.\d+", file)[0] ) )  
        with open(file, 'r') as f:
            for line in f:
                if (re.match('mean', line)): 
                    Rg.append( float(re.findall("\d+\.\d+", line)[0]) )
                    break 
            f.close() 

# print(T)
# print(Cv)

T = np.asarray(T)
Rg = np.asarray(Rg) 

Rg = [x for _, x in sorted(zip(T, Rg))] 
T = sorted(T) 

plt.figure()
plt.xlabel("Temperature (in E/$k_b$)")          
plt.ylabel("Radius of gyration")
plt.plot(T, Rg, linewidth=1, marker='^')

plt.show() 
plt.savefig("radius_of_gyration.png", dpi=1200) 

g = open('temp_rg.txt', 'w') 
g.write('T | Rg\n')
for i in range(len(T)): 
    g.write("{:0.2f} | {:0.2f}\n".format(T[i], Rg[i]) ) 
g.close() 
        
