#!/Users/satyend/opt/anaconda3/envs/CG/bin/python


import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
	col1 = "darkred" # "#2E93B7"
	col2 = "white" # "#B7532E"
	col3 = "darkblue" # "#B7532E"
	# col3 = "#807F7F"

	fig = plt.figure()
	figlegend = plt.figure(figsize=(3.25, 3))
	ax = fig.add_subplot(111)
	lines = []
	lines.append(ax.plot(range(10), np.random.randn(10), color=col1, markeredgecolor='k', ls='--', marker='o', markersize=8)[0])
	lines.append(ax.plot(range(10), np.random.randn(10), color=col2, markeredgecolor='k', ls='--', marker='o', markersize=8)[0])
	lines.append(ax.plot(range(10), np.random.randn(10), color=col3, markeredgecolor='k', ls='--', marker='o', markersize=8)[0])

	figlegend.legend(
		lines,
		[r'T=0.9', r'T=1.0',  r'T=1.1'],
		loc='center',
		frameon=False,
		ncol=1
	)
	figlegend.savefig('legend.png', bbox_inches='tight', dpi=1200)

