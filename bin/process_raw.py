#!/usr/bin/env python
import sys
import pandas as pd

### HUMAN PROTEIN ATLAS
## PART1: Opens the file and take the list of genes and tissues.
expression = dict()
tissues = set()
with open("../Reference/HPA_rna_tissue.tsv", 'r') as tis:
    tis = tis.read().splitlines()
    for lin in tis[1:100]:
        lin = lin.split('\t')
        tissues.add(lin[2])
    tissues = sorted(list(tissues))

    for line in tis[1:]:
        line = line.split('\t')
        key = line[1]
        if not key in expression:
            expression[key] = {}
            for z in tissues:
                tempdict = dict()
                tempdict[z] = 0
                expression[key].update(tempdict)

        tissu, val = line[2], line[3]
        if tissu in expression[key]:
            minidict = dict()
            minidict[tissu] = float(val)
            expression[key].update(minidict)

## PART2: Creates a pandas dataframe from the nested dictionary
## and write into an output file
df = pd.DataFrame.from_dict(expression,orient='index')
df.to_csv(path_or_buf='../Data/HPA_processed.tsv', sep='\t')

### GTEx
## Much simpler because already in a table-like format.
proc_gtex = open("../Data/GTEx_processed.tsv", 'w')
with open("../Reference/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", 'r') as gtex:
    gtex = gtex.read().splitlines()
    header = gtex[2].replace("gene_id\t", "").replace("Description\t", "")
    header = header.split('\t')
    proc_gtex.write('\t'.join(header)+'\n')
    for line in gtex[3:]:
        line = (line.split('\t'))[1:]
        proc_gtex.write('\t'.join(line) + '\n')


print ("Done!")
