#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
if __name__=="__main__":

    hexcolor_cg = 'steelblue'
    hexcolor_cc = 'lightskyblue'
    hexcolor_gg = 'lightcoral'
    hexcolor_gc = 'darkred'

    fig = pylab.figure()
    figlegend = pylab.figure(figsize=(2,1.5))
    ax = fig.add_subplot(111)
    lines = []
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_cc, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_cg, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gg, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gc, edgecolor='k'))

    figlegend.legend(lines, ('$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\parallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, s)$', '$(\\gamma,\\alpha _j) \\rightarrow (\\parallel, s)$'), 'center', frameon=False, ncol=4)
    figlegend.savefig('legend.png', dpi=2000, bbox_inches="tight")
