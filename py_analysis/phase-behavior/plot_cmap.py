import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mpl 
import numpy as np

def map_maker (hex_code):

    rgbcolor    = colors.hex2color (hex_code)
    cmap_colors = [(1,1,1), rgbcolor]
    my_cmap     = colors.LinearSegmentedColormap.from_list ('custom_cmap', cmap_colors)

    return my_cmap 


color_list = ['#1FB967', '#369DE8', '#B91F72', '#B9B41F']

for col in color_list:
    fig = plt.figure()
    ax  = plt.axes()
    
    my_cmap = map_maker (col)
    sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=mpl.colors.Normalize (vmin=0, vmax=1) )

    cbar = plt.colorbar (sm, orientation='vertical', format='%1.1f', shrink=0.5)
    cbar.set_ticks ( [0, 1] )
    
    cbar.set_ticklabels ( ["0.1", "6000"] )
    cbar.ax.tick_params (labelsize=10)
    for tick in cbar.ax.yaxis.get_major_ticks():
        tick.label2.set_fontweight('bold')
    
    ax.remove()

    plt.savefig (f"cmap_{color_list.index(col)}", bbox_inches="tight", dpi=1200)