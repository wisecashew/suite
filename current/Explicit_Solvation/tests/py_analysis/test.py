import numpy as np
import tidynamics as ti 

v = [[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]]

v_acf = ti.acf(v) 

print (v_acf)
print (len(v_acf))


