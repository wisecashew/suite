import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib 
from matplotlib.ticker import MultipleLocator
matplotlib.use("agg")
import argparse 
import time

parser = argparse.ArgumentParser(description="Plot the binodals from the .dump files.")
parser.add_argument("--dump",  dest="dump", nargs='+', type=str, action='store', help="Names of dump files.")
parser.add_argument("--csv",   dest="csv",  nargs='+', type=str, action='store', help="Names of csv files.")
parser.add_argument ("--draw-ticklabels",  dest='draw_tl',   action='store_true', help="Enter option to draw the ticklabels.", default=False)
parser.add_argument("--label", dest="label", type=str, action='store', help="Label of binodal.")
parser.add_argument("--xlim",  dest='xlim',  type=float, nargs='+', action='store', help="Provide a limits for x-axis.", default=[0,1])
parser.add_argument("--ylim",  dest='ylim',  type=float, nargs='+', action='store', help="Provide a limits for y-axis.", default=[0,1])
parser.add_argument ("--img",   dest="img",   type=str,   action='store', help="Name of image.", default="bintest.png")
args = parser.parse_args()

if __name__=="__main__":

    if len(args.dump) != len(args.csv):
        print(f"Number of dump files and number of csv files are not congruent.", flush=True)
        exit()

    start = time.time()
    if args.label == "UCST":
        colors = ["springgreen", "mediumseagreen", "seagreen"]
    elif args.label == "LOOP":
        colors = ["gold", "goldenrod", "darkgoldenrod"]
    elif args.label == "NECK":
        colors = ["coral", "red", "firebrick"]
    markers = ["^", "o", "s"]

    fig = plt.figure(num=0, figsize=(4,3), )
    ax  = plt.axes()
    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True, which='both')

    for i in range(len(args.csv)):
        try:
            df = pd.read_csv(args.csv[i], sep=',', engine="python", names=["phi", "T"])
            Texp_max = np.max(df["T"].values)
            ax.scatter(df["phi"].values, df["T"].values, marker=markers[i%len(markers)], s=50, c=colors[i], edgecolors='k', label="experimental data", zorder=25)
        except:
            print(f"No csv file.")

    for i in range(len(args.dump)):
        df = pd.read_csv(args.dump[i], engine='python', sep='\|', names=["left", "right", "T"], header=0)
        a_left = df["left"].values 
        a_right = df["right"].values
        a_temp  = df["T"].values

        temp_mask = np.logical_and(a_temp>=args.ylim[0], a_temp<=args.ylim[1])
        a_left    = a_left[temp_mask]
        a_right   = a_right[temp_mask]
        a_temp    = a_temp[temp_mask]

        left_mask = np.logical_and(a_left>=args.xlim[0], a_left<=args.xlim[1])
        a_left    = a_left[left_mask]
        a_right   = a_right[left_mask]
        a_temp    = a_temp[left_mask]

        right_mask = np.logical_and(a_right>=args.xlim[0], a_right<=args.xlim[1])
        a_left    = a_left[right_mask]
        a_right   = a_right[right_mask]
        a_temp    = a_temp[right_mask]

        print(f"left  = {a_left}")
        print(f"right = {a_right}")
        print(f"temp  = {a_temp}")

        ax.plot(a_left,  a_temp, c=colors[i], mec='k', lw=3)
        ax.plot(a_right, a_temp, c=colors[i], mec='k', lw=3)

    ax.set_ylim(args.ylim[0], args.ylim[1])
    ax.set_xlim(args.xlim[0], args.xlim[1])
    ax.set_yticks(np.linspace(args.ylim[0], args.ylim[1], 5))
    ax.set_xticks(np.linspace(args.xlim[0], args.xlim[1], 5))
    # ax.minorticks_on()
    # Set the number of minor ticks
    minor_ticks_between_major = 4
    nticks = 5

    # Set minor tick locators
    # ax.xaxis.set_minor_locator(MultipleLocator(5))
    # ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator((args.xlim[1] - args.xlim[0]) / (nticks - 1) / (minor_ticks_between_major + 1)))
    ax.yaxis.set_minor_locator(plt.MultipleLocator((args.ylim[1] - args.ylim[0]) / (nticks - 1) / (minor_ticks_between_major + 1)))

    if args.draw_tl:
        pass
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    fig.savefig(args.img, dpi=1200, bbox_inches="tight")
    stop = time.time()
    print(f"Time for computation is {stop - start} seconds.", flush=True)

