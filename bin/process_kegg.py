#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys
"""
Creates a ENSG to Gene Symbol dictionary from Ensembl Database (downloaded
10 July 2018). Then converts KEGG dB to gene symbol for further processing.
"""

print ("Preparing ENSG conversion dictionary...")
ensg = dict()
with open("../Reference/ensembl_ensg.tsv", 'r') as tis:
    for line in tis:
        line = line.split('\t')
        if not line[0] in ensg:
            ensg[line[0]] = line[1].rstrip('\n')

print ("Converting KEGG database to gene symbols...")
cluster = dict()
with open("../Reference/170508_kegg_short.tsv", 'r') as kegg:
    no_gs = open("../Reference/NO_ENSG_GeneSymb.tsv", 'w')
    for line in kegg:
        line = line.replace("_-_HOMO_SAPIENS_(HUMAN)\n", "")
        line = line.split('\t')
        if line[0] not in ensg:
            no_gs.write('\t'.join(line)+'\n')
        else:
            line[0] = ensg[line[0]]
            if not line[1] in cluster:
                cluster[line[1]] = set([line[0]])
            else:
                cluster[line[1]].update(set([line[0]]))

with open("../Reference/KEGG_GS.tsv", 'w') as kgs:
    for clust in cluster:
        values = cluster[clust]
        kgs.write(clust+'\t')
        kgs.write('\t'.join(values)+'\n')


print ("Done!")
