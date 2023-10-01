#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
if __name__=="__main__":

    fig = pylab.figure()
    figlegend = pylab.figure(figsize=(50,25))
    ax = fig.add_subplot(111)
    lines = []
    lines.append(ax.scatter(range(10), pylab.randn(10), color="coral", edgecolor='k'))
    lines.append(ax.scatter(range(10), pylab.randn(10), color="white", edgecolor='k'))
    # lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_cg, edgecolor='k'))
    # lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gg, edgecolor='k'))
    # lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gc, edgecolor='k'))

    figlegend.legend(lines, ('$\\exists \\varphi _{p}, \\varphi _{s}, \\varphi _{c}: D<0$', "stable"), frameon=False, ncol=1) #  ('$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\parallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, s)$', '$(\\gamma,\\alpha _j) \\rightarrow (\\parallel, s)$'), 'center', frameon=False)
    figlegend.savefig('legend.png', bbox_inches="tight", dpi=1200)
