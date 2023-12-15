#!/Users/satyend/opt/anaconda3/envs/CG/bin/python

import pylab
if __name__=="__main__":

    hexcolor_cg = '#369DE8' # 'steelblue'
    hexcolor_cc = '#1FB967' # 'lightskyblue'
    hexcolor_gg = '#B9B41F' # 'lightcoral'
    hexcolor_gc = '#B91F72' # 'darkred'

    fig = pylab.figure()
    figlegend = pylab.figure(figsize=(50,25))
    ax = fig.add_subplot(111)
    lines = []
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_cg, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_cc, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gg, edgecolor='k'))
    lines.append(ax.bar(range(10), pylab.randn(10), color=hexcolor_gc, edgecolor='k'))

    figlegend.legend(lines, ('$\\mathcal{R}_0$', '$\\mathcal{R}_1\'$', '$\\mathcal{R}_2$', '$\\mathcal{R}_3\'$'), 'center', frameon=False, ncol=4)
    figlegend.savefig('legend.png', bbox_inches="tight", dpi=1200)
