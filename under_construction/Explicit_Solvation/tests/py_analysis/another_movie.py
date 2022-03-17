import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
from matplotlib.animation import FuncAnimation

# setting up some subplots
fig, ax = plt.subplots(figsize=(5,3)) 
ax.set(xlim=(-3,3), ylim=(-1,1))

x = np.linspace(-3,3,91)
t = np.logspace(4,100, 91)

F = np.sin(x+t)

scat = ax.scatter(x, F)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'))

def animate(i):
	y = F + np.cos(i) 
	scat = ax.scatter(x,y)

# FuncAnimation requires two arguments: the figure object ‘fig’ and 
# animation function ‘animate’. The interval argument is optional and 
# sets the interval between frames in milliseconds. The frames 
# argument is needed only if the animation is to be exported. 
# The final two lines may or may not be necessary depending on whether 
# you’re working interactively (e.g. with IPython).
ani = FuncAnimation(scat, animate, interval=100, frames=16)
ani.save('movie2.mp4', writer=writer)