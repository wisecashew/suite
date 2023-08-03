import matplotlib.pyplot as plt
import matplotlib as mpl 
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mplc
from matplotlib.colors import LinearSegmentedColormap

# Define the color gradients and their corresponding positions
c1 = "#1AC6E5"
c2 = "#244BDB"
c3 = "#C157DB"
c4 = "#D95B59"
c5 = "#BFDA58"
c6 = "#2ED16A"
colors = [c1, c2, c2, c3, c3, c4, c4, c5, c5, c6] 
positions_raw   = np.array([-1, -0.53, -0.53, -0.525, -0.5250, -0.5235, -0.5235, -0.5, -0.5, 0]) 
k = 1/(np.max(positions_raw) - np.min (positions_raw))
d = 0-np.min(positions_raw)
positions_clean = k * ( positions_raw + d ) 

# Create the colormap
cmap = LinearSegmentedColormap.from_list('my_colormap', list(zip(positions_clean, colors)))

norm = mplc.Normalize (vmin=-1, vmax=0)

plt.rcParams['font.family'] = 'Arial'
font = {'color':  'black',
        'weight': 'normal'}

norm = mpl.colors.Normalize ( vmin=-1, vmax=0)

fig, ax = plt.subplots (figsize=(1.0,5.0)) # (figsize=(1,4))

# my_cmap = map_maker (col)
sm = plt.cm.ScalarMappable ( cmap=cmap, norm=norm )

cbar = plt.colorbar (sm, orientation='vertical', aspect=25, format='%1.1f', pad=0.15)
# cbar.ax.set_size_inches (2.5,1)
cbar.set_ticks ( [-1, -0.5, 0] )

# cbar.set_ticklabels ( ["-10", "-3", "-1", "-0.3", "-0.1"] )
cbar.ax.tick_params (labelsize=8)

ax.remove()

plt.savefig (f"cmap_gradient", bbox_inches="tight", dpi=1200)
