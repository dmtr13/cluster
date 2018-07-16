#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
    print ("Processing metric: {}".format(metrics))
    files = glob.glob("../Data/out*_{}*I20".format(metrics))
    files.sort()
    aggregate = [] #[threshold, [binnings]]
    for th, f in enumerate(files):
        size = []
        with open(f, 'r') as fil:
            for line in fil:
                line = line.rstrip('\n').split('\t')
                size.append(len(line))
        aggregate.append([threshold[th], binning(size)])
    return aggregate

def plot(aggregate, metrics):
    print ("Plotting...")

    N = len(bin)
    ind = np.arange(N)
    fig, ax = plt.subplots()
    width = 1/16
    p0 = ax.bar(ind-4*width, aggregate[0][1], width, color='b')
    p1 = ax.bar(ind-3*width, aggregate[1][1], width, color='g')
    p2 = ax.bar(ind-2*width, aggregate[2][1], width, color='r')
    p3 = ax.bar(ind-width, aggregate[3][1], width, color='c')
    p4 = ax.bar(ind, aggregate[4][1], width, color='m')
    p5 = ax.bar(ind+width, aggregate[5][1], width, color='k')
    p6 = ax.bar(ind+2*width, aggregate[6][1], width, color='orange')
    p7 = ax.bar(ind+3*width, aggregate[7][1], width, color='steelblue')
    p8 = ax.bar(ind+4*width, aggregate[8][1], width, color='brown')

    title = 'Cluster_Size_Distribution_{}'.format(metrics)
    ax.set_title('Distribution of Cluster Size for {} Metric'.format(metrics))
    ax.legend((threshold[1:6]))
    ax.set_xticks(ind+width/N)
    ax.set_xticklabels(bin, rotation=25)
    ax.set_ylabel('Count')
    ax.set_xlabel('Cluster Size')
    ax.autoscale_view()
    plt.savefig("../Data/{}.pdf".format(metrics))
    return True

filepath = "../Data/MCL_Cluster_Distribution.xlsx"
writer = pd.ExcelWriter(filepath, engine='xlsxwriter')

from collections import defaultdict
pt_dict = defaultdict(list)

for z, i in enumerate(metrics):
    aggregate = clust_th(i)
    plot(aggregate, i)

    ### per threshold ##########################################################
    for th in aggregate:
        if th[0] in pt_dict:
            pt_dict[th[0]].append([th[1]])
        else:
            pt_dict[th[0]] = [th[1]]
    ############################################################################

    ### per metric #############################################################
    if z == len(metrics)-1:
        print ("\nAggregating cluster sizes per metric...")
    pm = np.empty((len(threshold), len(bin)), dtype=int)
    for row in range(len(threshold)):
        pm[row, ] = np.array(aggregate[row][1])
    pm = pd.DataFrame(pm, index=threshold, columns=bin)
    pm.to_excel(writer, sheet_name=i)
    ############################################################################

print ("Aggregating cluster sizes per threshold...")
for th in pt_dict:
    thres = np.empty((len(metrics), len(bin)), dtype=int)
    for row in range(len(metrics)):
        thres[row, ] = np.array(pt_dict[th][row])
    thres = pd.DataFrame(thres, index=metrics, columns=bin)
    thres.to_excel(writer, sheet_name=str(th))

writer.close()
print ("Done!\n")
