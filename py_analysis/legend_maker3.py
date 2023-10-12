import pylab

if __name__ == "__main__":
    fig = pylab.figure()
    figlegend = pylab.figure(figsize=(50, 25))
    ax = fig.add_subplot(111)
    lines = []
    lines.append(ax.scatter(range(10), pylab.randn(10), color="silver"))
    lines.append(ax.scatter(range(10), pylab.randn(10), color="darkred"))
    lines.append(ax.scatter(range(10), pylab.randn(10), color="blue", marker='o'))
    lines.append(ax.scatter(range(10), pylab.randn(10), color="red", marker='o'))
    plot_line = ax.plot(range(10), pylab.randn(10), color='skyblue', lw=1, linestyle='--')
    lines.append(plot_line)

    legend_labels = ('binodal arm', 'conjugate binodal arm', 'unstable', 'stable', 'tieline')

    # Create the legend for scatter plots and add it to figlegend
    figlegend.legend(lines[:4], legend_labels[:4], frameon=False, ncol=1)

    # Create the legend for the line plot (ax.plot) and add it to figlegend
    figlegend.legend([plot_line[0]], [legend_labels[4]], frameon=False, loc='upper right')

    # Adjust the position of the legends
    figlegend.subplots_adjust(top=0.9, right=0.8)

    # Save or display the legends
    figlegend.savefig('legends.png', bbox_inches='tight', transparent=True)

