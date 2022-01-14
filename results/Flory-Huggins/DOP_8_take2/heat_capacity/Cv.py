import os 
import re 
import matplotlib.pyplot as plt
import numpy as np 

cwd = os.getcwd() 
files = os.listdir() 

T = []
Cv =[] 

for file in files:
    if (re.match('heat_cap', file)):
        # print(re.findall("\d+\.\d+", file))
        T.append( float(re.findall("\d+\.\d+", file)[0] ) )  
        with open(file, 'r') as f:
            for line in f:
                if (re.match('heat capacity', line)): 
                    Cv.append( float(re.findall("\d+\.\d+", line)[0]) )
                    break 
            f.close() 

# print(T)
# print(Cv)

T = np.asarray(T)
Cv = np.asarray(Cv) 

Cv = [x for _, x in sorted(zip(T, Cv))] 
T = sorted(T) 

plt.figure()
plt.xlabel("Temperature (in $k_b$T)")          
plt.ylabel("Heat capacity (in $k_b$)")
plt.plot(T, Cv, linewidth=1, marker='^')

plt.savefig("transition.png", dpi=1200) 
plt.show() 

g = open('temp_cv.txt', 'w') 
g.write('T | Cv\n')
for i in range(len(T)): 
    g.write("{:0.2f} | {:0.2f}\n".format(T[i], Cv[i]) ) 
g.close() 
        
