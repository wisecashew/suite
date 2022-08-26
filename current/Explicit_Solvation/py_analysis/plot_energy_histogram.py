#!/usr/licensed/anaconda3/2020.7/bin/python

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt 
import sys 

parser = argparse.ArgumentParser (description="Read the energydump, and plot the energy histogram") 
parser.add_argument ('-i', nargs='+', metavar='energydump', dest='i', type=str, action='store', help="Address to energy dump file to parse information.", default="energydump") 
parser.add_argument ('-s', metavar='S', dest='s', type=int, action='store', help='Starting index to construct histogram.', default=0)
parser.add_argument ('-T', metavar='T', dest='T', type=float, action='store', help='Input temperature of system') 

args = parser.parse_args() 

if __name__ == "__main__": 
    
    fig = plt.figure()
    ax  = plt.axes() 
    labels = ["biased", "unbiased"]
    j = 1
    colors = ["red", "darkblue"] 
    edict = {} 
    k = 1
    for (legend, address, c) in zip(labels, args.i, colors):
        df = pd.read_csv(address, sep=' \| ', names=["energy", "mm_tot", "mm_aligned", "mm_naligned", "ms_tot", "ms_aligned", "ms_naligned", "time_step"], engine="python", skiprows=args.s)
    
        energy_store = df["energy"].values 
        h = plt.hist (energy_store, bins=20, density=True, alpha=j, label=legend, ls='dashed', color=c, edgecolor='k')
        j = j/2 
        edict[k] = energy_store 
        k += 1

    
    plt.legend(labels)
    plt.savefig('energy_histogram_'+str(args.T)+'.png', dpi=1200) 
    plt.show() 

    print ("Moving on to calculate KL divergence...")

    # partition the energy store into a bunch of bins 
    # p2 is the true distribution
    true_min = np.min([np.min(edict[1]), np.min(edict[2])])
    true_max = np.max([np.max(edict[1]), np.max(edict[2])])
    rnge = (true_max-true_min) 
    nbins = 20 
    bwidth = 1/nbins 
    p1 = np.zeros(nbins) 
    p2 = np.zeros(nbins)
    
    for elem in edict[1]:
        ratio = (elem-true_min)/rnge 
        idx   = ratio//bwidth 
        # print (idx)
        # print (int(idx))
        idx = int(idx)
        if idx == nbins:
            idx = nbins-1
        p1[idx] += 1


    for elem in edict[2]:
        ratio = (elem-true_min)/rnge 
        idx   = ratio//bwidth 
        #  print (idx)
        # print (int(idx))
        idx = int(idx)
        if idx == nbins:
            idx = nbins-1
        p2[idx] += 1
    
    print (p1)
    print (p2)
    p1 = p1/np.sum(p1)
    p2 = p2/np.sum(p2)
    idx_to_del = [] 
    for i in range (len(p1)):
        if ( p1[i] == 0 or p2[i] == 0 ):
            if (abs(p1[i]-p2[i]) > 10 ):
                print ("Difference too large. Maybe a problem.")
                print ("Index is ", i)
                sys.exit("Issue detected.") 
            idx_to_del.append(i) 

    s = 0
    for x in idx_to_del:
        p1 = np.delete(p1, x-s)
        p2 = np.delete(p2, x-s)
        s = s+1

    kl_div = np.sum( p1*np.log(p1/p2) ) 

    print ("KL divergence = ", kl_div)

    
