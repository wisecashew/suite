import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl 
import numpy as np

fig = plt.figure()
ax  = plt.axes()
my_cmap = cm.PiYG
sm = plt.cm.ScalarMappable ( cmap=my_cmap, norm=mpl.colors.TwoSlopeNorm (vmin=-10, vcenter=-5, vmax=-0) )

cbar = plt.colorbar (sm, orientation='vertical', format='%1.1f', shrink=0.5)
cbar.set_ticks ( [-10, -7.5, -5, -2.5, 0] )
# cbar.set_ticklabels ( ["$-\mathbf{10^{-1}}$", "-$\mathbf{10^{-2}}$","-$\mathbf{10^{-3}}$",  "$\mathbf{0}$", "$\mathbf{10^{-3}}$", "$\mathbf{10^{-2}}$", "$\mathbf{10^{-1}}$" ] )
cbar.set_ticklabels ( ["-10", "-5", "-1", "-0.75", "-0.1"] )
cbar.ax.tick_params (labelsize=10)
for tick in cbar.ax.yaxis.get_major_ticks():
    tick.label2.set_fontweight('bold')
# cbar.ax.ticklabel_format (labelweight='bold')
ax.remove()

plt.savefig ("cmap_piyg", bbox_inches="tight", dpi=1200)