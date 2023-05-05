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


plt.rcParams['font.family'] = 'Arial'
font = {'color':  'black',
        'weight': 'normal'}
color_list = ['#1FB967', '#369DE8', '#B91F72', '#B9B41F']

emsa_list = [-10, -5, -1, -0.5, -0.1]
emsa_list = [-10, -5, -1, -0.5, -0.1] # [-0.6, -0.4] # np.arange (elow, ehigh+0.05, 0.05)
elow    = np.min(emsa_list)
ehigh   = np.max (emsa_list)
ecenter = (elow+ehigh)/2
norm = mpl.colors.Normalize ( vmin=-4, vmax=0)


for col in color_list:
    fig, ax = plt.subplots (figsize=(1.0,5.0)) # (figsize=(1,4))

    # my_cmap = map_maker (col)
    sm = plt.cm.ScalarMappable ( cmap=cm.turbo, norm=norm )

    cbar = plt.colorbar (sm, orientation='vertical', aspect=25, format='%1.1f', pad=0.15)
    # cbar.ax.set_size_inches (2.5,1)
    cbar.set_ticks ( [-4, -2, -0] )
    
    # cbar.set_ticklabels ( ["-10", "-3", "-1", "-0.3", "-0.1"] )
    cbar.ax.tick_params (labelsize=8)

    # for tick in cbar.ax.yaxis.get_major_ticks():
        # tick.label2.set_fontweight('bold')
    
    ax.remove()

    plt.savefig (f"cmap_turbo", bbox_inches="tight", dpi=1200)
