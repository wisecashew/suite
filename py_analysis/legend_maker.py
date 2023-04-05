#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
fig = pylab.figure()
figlegend = pylab.figure(figsize=(3,2))
ax = fig.add_subplot(111)
lines = []
lines.append(ax.bar(range(10), pylab.randn(10), color='darkred', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='lightcoral', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='steelblue', edgecolor='k'))
lines.append(ax.bar(range(10), pylab.randn(10), color='lightskyblue', edgecolor='k'))

figlegend.legend(lines, ('$\\mathbf{i \\rightarrow m ^{\\parallel} }$', '$\\mathbf{ i \\rightarrow m }^{\\perp} $', '$\\mathbf{ i \\rightarrow s ^ {\\parallel} }$', '$\\mathbf{ i \\rightarrow s }^{\\perp} $'), 'center', frameon=False)
# fig.show()
# figlegend.show()
figlegend.savefig('legend.png', dpi=1200, bbox_inches="tight")