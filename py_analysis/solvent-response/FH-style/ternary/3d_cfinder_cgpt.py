import numpy as np

def f(x, y, c):
    z = np.sin(x+y) * c[:, :, np.newaxis]
    return z

x = np.arange(1,10)
y = np.arange(11,20)

c = np.arange(1,5).reshape((2,2))

print (f"c.shape = {c.shape}")
print (f"x.shape = {x.shape}")
print (f"y.shape = {y.shape}")
z = f(x, y, c)
print (f"z.shape = {z.shape}")
print (z)