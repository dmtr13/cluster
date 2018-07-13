#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import glob, math, subprocess, sys
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

"""
Non-modular script to calculate the statistics of the clusters as
generated by MCL.
"""

threshold = [0.9, 0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999]
metrics = ["Relative-Euclidean", "Euclidean",
            "Mass-Distance", "Manhattan"]

def cluster(glob_pattern):
    files = glob.glob("../Data/out*_{}*I20".format(glob_pattern))
    files.sort()
    assert len(files) == len(threshold), "Unequal #Threshold & #Files"
    files.sort()
    aggregate = []
    for th, f in enumerate(files):
        n_cluster = 0
        cluster_size = []
        with open(f, 'r') as cl:
            cl = cl.readlines()
            for ea_cl in cl:
                n_cluster += 1
                cluster_size.append(len(ea_cl.split('\t')))
        assert n_cluster == len(cluster_size), "#Cluster != #Cluster_size"
        counter = Counter()
        for gene in cluster_size:
            counter[gene] += 1
        aggregate.append([n_cluster, cluster_size, counter, th])
    return aggregate

def clus_no(aggr):
    thelist = [item[0] for item in aggr]
    # return thelist[1:4]
    return thelist

def clus_size(aggr):
    flatlist = dict()
    for thresh in aggr:
        th = thresh[3]
        counter = [(k,v) for k,v in thresh[2].items()]
        for counts in counter:
            clsize, freq = int(counts[0]), int(counts[1])
            flatlist.append([float(threshold[th]), clsize, freq])
    return flatlist

print ("Processing statistics...")
rel_euc = cluster(metrics[0])
euc = cluster(metrics[1])
md = cluster(metrics[2])
man = cluster(metrics[3])

### No. of Clusters vs. Threshold ####################################
print ("Plotting distribution of no. of clusters against cut-off...")
##
# threshold = threshold[1:4] ##Also check f(clus_no)!
##
N = len(threshold)
ind = np.arange(N)
fig, ax = plt.subplots()
width = 0.225
p1 = ax.bar(ind-width, clus_no(rel_euc), width, color='b')
p2 = ax.bar(ind, clus_no(euc), width, color='r')
p3 = ax.bar(ind+width, clus_no(md), width, color='g')
p4 = ax.bar(ind+2*width, clus_no(man), width, color='k')

ax.set_title('Distribution of No. of Clusters Across Different Thresholds')
ax.set_xticks(ind + width / N)
ax.set_xticklabels(threshold)

ax.legend((p1[0], p2[0], p3[0], p4[0]), metrics)
ax.set_ylabel("No. of Clusters")
ax.set_xlabel("Threshold")
ax.autoscale_view()
plt.savefig("../Data/Threshold_vs_NCluster.pdf")
plt.close()
######################################################################

metr = {"Relative-Euclidean":clus_size(rel_euc),
            "Euclidean":clus_size(euc),
            "Mass-Distance":clus_size(md),
            "Manhattan":clus_size(man)}

# print ("Processing distribution of no. of genes per cluster...")
# for fi in metrics:
#     print ("Metric:\t{}".format(fi))
#     with open("../Data/{}.stats".format(fi), 'w') as st:
#         st.write("Metric:\t{}\n".format(fi))
#         st.write("Threshold | Cluster Size (genes)| Count\n")
#         data = metr[fi]
#         for line in data:
#             st.write(' | '.join(line)+"\n")

print ("Done!\n")
