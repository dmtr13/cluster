#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, glob
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

metrics = ["Relative-Euclidean", "Euclidean",
            "Mass-Distance", "Manhattan"]
threshold = [0.9, 0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999]
def clust_th(metrics):
    print ("Metric: {}".format(metrics))
    files = glob.glob("../Data/out*_{}*I20".format(metrics))
    files.sort()
    # print (files)
    aggregate = [] #[threshold, [binnings]]
    for th, f in enumerate(files):
        size = []
        with open(f, 'r') as fil:
            for line in fil:
                line = line.rstrip('\n').split('\t')
                size.append(len(line))
        aggregate.append([threshold[th], binning(size)])
    # print (aggregate)

    print ("Plotting dist. of cluster size across different thresholds...")

    N = len(bin)
    ind = np.arange(N)
    fig, ax = plt.subplots()
    width = 1/16
    # p0 = ax.bar(ind-4*width, aggregate[0][1], width, color='b')
    p1 = ax.bar(ind-3*width, aggregate[1][1], width, color='g')
    p2 = ax.bar(ind-2*width, aggregate[2][1], width, color='r')
    p3 = ax.bar(ind-width, aggregate[3][1], width, color='c')
    p4 = ax.bar(ind, aggregate[4][1], width, color='m')
    p5 = ax.bar(ind+width, aggregate[5][1], width, color='k')
    # p6 = ax.bar(ind+2*width, aggregate[6][1], width, color='orange')
    # p7 = ax.bar(ind+3*width, aggregate[7][1], width, color='steelblue')
    # p8 = ax.bar(ind+4*width, aggregate[8][1], width, color='brown')

    title = 'Cluster_Size_Distribution_{}'.format(metrics)
    ax.set_title('Distribution of Cluster Size for {} Metric'.format(metrics))
    ax.legend((threshold[1:6]))
    ax.set_xticks(ind+width/N)
    ax.set_xticklabels(bin, rotation=25)
    ax.set_ylabel('Count')
    ax.set_xlabel('Cluster Size')
    ax.autoscale_view()
    plt.savefig("../Data/{}.pdf".format(metrics))
    # plt.show()

    return True

for i in metrics:
    clust_th(i)
# clust_th("Euclidean")
print ("Done")

# def clus_size(pathtofile):
#     ### OUTFILE ################################################################
#     filename_extension = os.path.basename(sys.argv[1])
#     directory = os.path.dirname(sys.argv[1])
#     fn, ext = os.path.splitext(filename_extension)
#     fn = fn.replace('out', '').replace('.tsv.I20', '')
#     fn = fn.split('_')
#     fn = ' '.join(fn[2:4])
#     # fn = fn.replace('out.HPA_processed_', '')
#     # fn = fn.replace('_MCL.tsv.I20', '')
#     print ("Processing {}...".format(fn))
#     # outfile = directory+'/'+fn+'_KEGG'+'.tsv'
#     ############################################################################
#
#     size = []
#     with open(pathtofile, 'r') as fil:
#         for line in fil:
#             line = line.rstrip('\n').split('\t')
#             size.append(len(line))
#     # print (size)
#     # plt.hist(size, bins=bin)
#     # plt.ylim([0, 30])
#     # plt.show()
#     y = binning(size)
#     x_shade = range(len(y))
#     plt.figure()
#     plt.bar(x_shade, y)
#     plt.xticks(x_shade, bin, rotation=20)
#     plt.ylabel('Counts')
#     plt.xlabel('Size of Cluster')
#     plt.tight_layout()
#     plt.show()
#     # size = [[x,size.count(x)] for x in set(size)]
#     # size = sorted(size, key=itemgetter(0))
#     # x, y = [], []
#     # for z in size:
#     #     x.append(str(z[0]))
#     #     y.append(z[1])
#     # x_shade = range(0, len(x)*4, 4)
#     # plt.figure(figsize=(6,3))
#     # plt.xticks(x_shade, x, rotation='vertical')
#     # plt.bar(x_shade, y)
#
#     # for y in x:
#     #     print (y)
#
# if __name__ == "__main__":
#     print(clus_size(sys.argv[1]))
