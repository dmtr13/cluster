#!/usr/bin/env python3
import sys, time
import pandas as pd

start = time.time()
"""
A non-modular script that converts the Human Protein Atlas' RNA counts from
37 different tissues and GTEx's RNA counts into a tabular format.
"""

### HUMAN PROTEIN ATLAS
print ("Processing HPA...")
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
print ("Processing GTEx...")
## Much simpler because already in a table-like format.
proc_gtex = open("../Data/GTEx_processed.tsv", 'w')
### GTEx has non-coding genes, so let's filter those
pcg = []
with open("../Reference/genesymb.txt", 'r') as gs:
    gs = gs.readlines()
    for line in gs[1:]:
        pcg.append(line.rstrip('\n'))
with open("../Reference/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", 'r') as gtex:
    gtex = gtex.read().splitlines()
    header = gtex[2].replace("gene_id\t", "").replace("Description\t", "")
    header = header.split('\t')
    proc_gtex.write('\t'.join(header)+'\n')
    temp = []
    npcg_count = 0
    for line in gtex[3:]:
        line = (line.split('\t'))[1:]
        if line[0] in pcg:
            temp.append('\t'.join(line))
        else:
            npcg_count += 1
    temp = sorted(temp)
    for i in temp:
        proc_gtex.write(i+'\n')
    print ("Excluded {} non-protein coding genes.".format(npcg_count))

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
