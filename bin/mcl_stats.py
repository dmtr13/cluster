#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import glob, math, subprocess, sys
import matplotlib.pyplot as plt
import numpy as np

## define a fx that takes the glob pattern for each metric and
## save the length of lines (no. of clusters) into a list
threshold = [0, 0.5, 0.9, 0.925, 0.95, 0.975, 0.99]
metrics = ["Relative-Euclidean", "Euclidean",
            "Mass-Distance", "Manhattan"]

def cluster(glob_pattern):
    files = glob.glob("../Data/2500/*_{}*_I20.tsv".format(glob_pattern))
    assert len(files) == len(threshold), "Unequal #Threshold & #Files"
    files.sort()
    aggregate = []
    for f in files:
        n_cluster = 0
        cluster_size = []
        with open(f, 'r') as cl:
            cl = cl.readlines()
            for ea_cl in cl:
                n_cluster += 1
                cluster_size.append(len(ea_cl.split('\t')))
        assert n_cluster == len(cluster_size), "#Cluster != #Cluster_size"
        aggregate.append([n_cluster, cluster_size])
    return aggregate

def clus_no(aggr):
    return [item[0] for item in aggr]

def clus_size(aggr):
    counts = []
    for item in aggr:
        counts.append([[x,item[1].count(x)] for x in set(item[1])])
    return counts

print ("Processing statistics...")
rel_euc = cluster(metrics[0])
euc = cluster(metrics[1])
md = cluster(metrics[2])
man = cluster(metrics[3])

### No. of Clusters vs. Threshold ####################################
print ("Plotting distribution of no. of clusters against cut-off...")
N = 7
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
plt.savefig("../Data/2500/Threshold_vs_NCluster.pdf")
plt.close()
######################################################################

metr = {"Relative-Euclidean":clus_size(rel_euc),
            "Euclidean":clus_size(euc),
            "Mass-Distance":clus_size(md),
            "Manhattan":clus_size(man)}

print ("Processing distribution of no. of genes per cluster...")
for fi in metrics:
    print ("Metric:\t{}".format(fi))
    with open("../Data/2500/{}.stats".format(fi), 'w') as st:
        st.write("Metric:\t{}\n".format(fi))
        st.write("Threshold\tCluster Size\tNo. of Genes\n")
        data = metr[fi]
        for th in range(len(data)):
            for size in data[th]:
                st.write(str(threshold[th])+'\t'+str(size[0]) \
                        +'\t'+str(size[1])+'\n')

print ("Done!\n")
