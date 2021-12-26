import os 
import re 
import matplotlib.pyplot as plt

cwd = os.getcwd() 
files = os.listdir() 

T = []
Cv =[] 

for file in files:
    if (re.match('heat_cap', file)):
        T.append( float(re.findall("\d+\.\d+", file)[0]) )  
        with open(file, 'r') as f:
            for line in f:
                if (re.match('heat capacity', line)): 
                    Cv.append( float(re.findall("\d+\.\d+", line)[0]) )
                    break 
            f.close() 

plt.xlabel("Temperature (in $k_b$T)")          
plt.ylabel("Heat capacity (in $k_b$)"
plt.scatter(T, Cv)

# plt.show() 
plt.savefig("heat_cap.png", dpi=1200) 

Cv = [x for _, x in sorted(zip(T, Cv))] 
T = sorted(T) 

g = open('cv.txt', 'w') 
g.write('T | Cv\n')
for i in range(len(T)): 
    g.write("{:0.2f} | {:0.2f}\n".format(T[i], Cv[i]) ) 
g.close() 
        
