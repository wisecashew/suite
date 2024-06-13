import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

# Generate random data for plots
x = np.linspace(0, 10, 100)
y1 = np.random.rand(100)
y2 = np.random.rand(100)
y3 = np.random.rand(100)

custom_font = FontProperties(fname="arial.ttf")

# Create plots
# plt.plot([], [], c="gold",      lw=2, label='tieline', ls='--')
# plt.plot([], [], c="orangered", lw=4, label="binodal arm (I)")
# plt.plot([], [], c="darkred",   lw=4, label="binodal arm (II)")
# plt.scatter([], [], c="red",       s=5, label="critical point")
# plt.scatter([], [], c="steelblue", s=5, label="metastable/stable")
# plt.scatter([], [], c="coral",     s=5, label="unstable")

plt.plot([], [], c='#8EDA25', ls='--', lw=1, markersize=2, marker='o', mec='k', label="solvent mixing")
plt.plot([], [], c='#258EDA', ls='--', lw=1, markersize=2, marker='o', mec='k', label="enthalpic bridging")
plt.plot([], [], c='#DA258E', ls='--', lw=1, markersize=2, marker='o', mec='k', label="geometric frustration")


# Display legend window only
plt.legend(loc='center', frameon=False, prop=custom_font)
plt.axis('off')
plt.savefig("legend.png", dpi=1200, bbox_inches='tight', transparent=True)

