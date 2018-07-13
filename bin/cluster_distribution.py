#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from operator import itemgetter

bin = ['1-10', '11-50', '51-100', '101-250', '250-500', '501-1000',
        '1001-2500', '>2500']
def binning(sizes):
    thesize = [0, 0, 0, 0, 0, 0, 0, 0]
    for sz in sizes:
        if sz <= 10:
            thesize[0] += 1
        elif sz <= 50:
            thesize[1] += 1
        elif sz <= 100:
            thesize[2] += 1
        elif sz <= 250:
            thesize[3] += 1
        elif sz <= 500:
            thesize[4] += 1
        elif sz <= 1000:
            thesize[5] += 1
        elif sz <= 2500:
            thesize[6] += 1
        else:
            thesize[7] += 1
    return thesize



def clus_size(pathtofile):
    ### OUTFILE ################################################################
    filename_extension = os.path.basename(sys.argv[1])
    directory = os.path.dirname(sys.argv[1])
    fn, ext = os.path.splitext(filename_extension)
    fn = fn.replace('out', '').replace('.tsv.I20', '')
    fn = fn.split('_')
    fn = ' '.join(fn[2:4])
    # fn = fn.replace('out.HPA_processed_', '')
    # fn = fn.replace('_MCL.tsv.I20', '')
    print ("Processing {}...".format(fn))
    # outfile = directory+'/'+fn+'_KEGG'+'.tsv'
    ############################################################################

    size = []
    with open(pathtofile, 'r') as fil:
        for line in fil:
            line = line.rstrip('\n').split('\t')
            size.append(len(line))
    # print (size)
    # plt.hist(size, bins=bin)
    # plt.ylim([0, 30])
    # plt.show()
    y = binning(size)
    x_shade = range(len(y))
    plt.figure()
    plt.bar(x_shade, y)
    plt.xticks(x_shade, bin, rotation=20)
    plt.ylabel('Counts')
    plt.xlabel('Size of Cluster')
    plt.tight_layout()
    plt.show()
    # size = [[x,size.count(x)] for x in set(size)]
    # size = sorted(size, key=itemgetter(0))
    # x, y = [], []
    # for z in size:
    #     x.append(str(z[0]))
    #     y.append(z[1])
    # x_shade = range(0, len(x)*4, 4)
    # plt.figure(figsize=(6,3))
    # plt.xticks(x_shade, x, rotation='vertical')
    # plt.bar(x_shade, y)

    # for y in x:
    #     print (y)

if __name__ == "__main__":
    print(clus_size(sys.argv[1]))
