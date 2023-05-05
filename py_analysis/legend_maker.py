#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
fig = pylab.figure()
figlegend = pylab.figure(figsize=(2,1.5))
ax = fig.add_subplot(111)
lines = []
lines.append(ax.bar(range(10), pylab.randn(10), color='darkred', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='lightcoral', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='steelblue', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='lightskyblue', edgecolor='k'))

figlegend.legend(lines, ('$ (\\gamma, \\alpha_j) \\rightarrow (\\parallel, s) $', '$ (\\gamma, \\alpha_j) \\rightarrow (\\nparallel, s) $', '$ (\\gamma, \\alpha_j) \\rightarrow (\\parallel, m) $', '$ (\\gamma, \\alpha_j) \\rightarrow (\\nparallel, m) $'), 'center', frameon=False, ncol=2)
# fig.show()
# figlegend.show()
figlegend.savefig('legend.png', dpi=1200, bbox_inches="tight")
