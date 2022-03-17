#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 23:06:40 2022

@author: satyendhamankar
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure( figsize=(10,6) ) 
ax = fig.add_subplot(111, projection='3d')

Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'))

def animate(i):
    x = np.linspace(0, 2*np.pi, 120)
    y = np.cos(i*x)
    z = x
    # ax.clear()
    print(i)
    im = ax.plot(x, y, z)

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=range(5), repeat=False, blit=False)
ani.save("movie.mp4", writer=writer)
plt.show()
