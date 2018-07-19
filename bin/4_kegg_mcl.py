#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os
import pandas as pd
import numpy as np

"""
Reads in KEGG dB with gene symbols instead of ENSGs and compare each pathway
to the clusters as determined by MCL.
"""

cluster = dict()
keys = []
with open("../Reference/KEGG_GS.tsv", 'r') as kgs:
    for line in kgs:
        line = line.replace("\n", '').split('\t')
        key = line[0]+("({})".format(len(line[1:])))
        keys.append(key)
        cluster[key] = [len(line[1:]), set(line[1:])]
# print (cluster)
### OUTFILE ####################################################################
filename_extension = os.path.basename(sys.argv[1])
directory = os.path.dirname(sys.argv[1])
fn, ext = os.path.splitext(filename_extension)
fn = fn.lstrip('out.').replace(".tsv", '')
print ("Processing {}...".format(fn))
outfile = directory+'/'+fn+'_KEGG'+'.tsv'
################################################################################

colcounts = []
mclgenes = dict()
inp = open(sys.argv[1], 'r')
for z, line in enumerate(inp):
    line = line.rstrip('\n').split('\t')
    if len(line) > 10:
        colcounts.append("Cluster#{a}-({b})".format(a=z+1, b=len(line)))
        mclgenes[z] = set(line)
array = np.zeros((len(keys), len(colcounts)), dtype=float)#, \
        #dtype=[('number', 'i4'), ('percentage', 'f6')])

kegg_counter = range(len(keys))
mcl_counter = range(len(colcounts))

def jaccard_sim(set1, set2):
    if len(set1) + len(set2) == 0:
        return 1.0
    else:
        index = len(cluster[keys[k]][1] & mclgenes[m]) /\
                len(cluster[keys[k]][1] | mclgenes[m])
        return index

for k in kegg_counter:
    for m in mcl_counter:
        # intersection = cluster[keys[k]][1] & mclgenes[m]
        # array[k, m] = intersection
        array[k, m] = jaccard_sim(cluster[keys[k]][1], mclgenes[m])


df = pd.DataFrame(array, index=keys, columns=colcounts)
df.to_csv(outfile, sep='\t')

print ("Done!")
