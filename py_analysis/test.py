import matplotlib.pyplot as plt
import pylab

if __name__ == "__main__":
    hexcolor_cg = 'steelblue'
    hexcolor_cc = 'lightskyblue'
    hexcolor_gg = 'lightcoral'
    hexcolor_gc = 'darkred'

    fig = plt.figure()
    figlegend = plt.figure(figsize=(2, 1.5))
    ax = fig.add_subplot(111)
    lines = []
    lines.extend(ax.plot(range(10), pylab.randn(10), color=hexcolor_cc))
    lines.extend(ax.plot(range(10), pylab.randn(10), color=hexcolor_cg))
    lines.extend(ax.plot(range(10), pylab.randn(10), color=hexcolor_gg))
    lines.extend(ax.plot(range(10), pylab.randn(10), color=hexcolor_gc))

    figlegend.legend(lines, ('$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\parallel, m)$', '$(\\gamma,\\alpha _j)\\rightarrow (\\nparallel, s)$', '$(\\gamma,\\alpha _j) \\rightarrow (\\parallel, s)$'), 'center', frameon=False, ncol=4)
    figlegend.savefig('legend.png', dpi=1200, bbox_inches="tight")

