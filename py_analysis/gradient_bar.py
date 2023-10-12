import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Define three colors
color3 = np.array([131, 159, 192]) / 255.0  
color2 = "white"
color1 = np.array([137, 245, 162]) / 255.0  

# Create a list of color points with positions and colors
colors = [(0.0, color1), (0.5, color2), (1.0, color3)]

# Create the custom colormap using LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)

# Create a sample gradient image using the colormap
gradient = plt.imshow([[0, 1]], cmap=custom_cmap)

# Add a colorbar to visualize the gradient colormap
plt.gca().set_visible(False)
cbar = plt.colorbar(gradient)
cbar.set_ticks ([0, 0.5, 1])
cbar.set_ticklabels([])

# Show the plot
plt.savefig("col_bar.png", bbox_inches="tight", dpi=1200)
