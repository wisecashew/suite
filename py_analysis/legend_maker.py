#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
if __name__=="__main__":

	col1 = "steelblue"
	col2 = "coral"
	col3 = "lavender"
	col4 = "seagreen"

	fig = pylab.figure()
	figlegend = pylab.figure(figsize=(50,25))
	ax = fig.add_subplot(111)
	lines = []
	lines.append(ax.plot(range(10), pylab.randn(10), color=col1, markeredgecolor='k', ls='--', marker='o')[0])
	lines.append(ax.plot(range(10), pylab.randn(10), color=col2, markeredgecolor='k', ls='--', marker='^')[0])
	lines.append(ax.plot(range(10), pylab.randn(10), color=col3, markeredgecolor='k', ls='--')[0])
	lines.append(ax.plot(range(10), pylab.randn(10), color=col4, markeredgecolor='k', ls='--', marker='o')[0])

	figlegend.legend(lines, ('$\\tilde{R}_g$', '$\\langle \\tilde{R}_g \\rangle_i$', '$\\langle \\tilde{R}_g \\rangle$', 'ACF'), 
					 loc='center', frameon=False, ncol=4)
	figlegend.savefig('legend.png', bbox_inches="tight", dpi=1200)

