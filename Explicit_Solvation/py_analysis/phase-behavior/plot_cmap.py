import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl 
import numpy as np

fig = plt.figure()
ax  = plt.axes()
my_cmap = cm.coolwarm
sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1) )

cbar = plt.colorbar (sm, orientation='vertical', format='%.0e', shrink=0.5)
cbar.set_ticks ( np.linspace(0, 1, 5))
# cbar.set_ticklabels ( ["$-\mathbf{10^{-1}}$", "-$\mathbf{10^{-2}}$","-$\mathbf{10^{-3}}$",  "$\mathbf{0}$", "$\mathbf{10^{-3}}$", "$\mathbf{10^{-2}}$", "$\mathbf{10^{-1}}$" ] )
cbar.set_ticklabels ( ["$\mathbf{10^{-2}}$", "$\mathbf{10^{-1}}$","$\mathbf{10^{0}}$",  "$\mathbf{10^1}$", "$\mathbf{10^{2}}$" ] )
cbar.ax.tick_params (labelsize=10)
for tick in cbar.ax.yaxis.get_major_ticks():
    tick.label2.set_fontweight('bold')
# cbar.ax.ticklabel_format (labelweight='bold')
ax.remove()

plt.savefig ("just_colormap", bbox_inches="tight", dpi=1200)