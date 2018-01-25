#!/bin/usr/env python3
import math, sys
import numpy as np


## this might be wrong, as we're looking into variations of the genes
## between the tissues, so the matrix should be between genes and genes
## instead of tissues and tissues


filename = open(sys.argv[1], 'r')
filename = filename.read().splitlines()
header = filename[0].split('\t')
del header[0]
c = len(header)
dimension = len(header)
genes = []
values = []

for lines in filename[1:]:
    lines = lines.split()
    genes.append(lines[0])
    values.append([float(x) for x in lines[1:]])
tissue_tpm = dict(zip(header, values))

def euclidean(tissue_tpm):
    eucvalues = []
    for t1 in tissue_tpm:
        print ("Calculating %s set..." % t1)
        for t2 in tissue_tpm:
            theset = [(a-b)**2 for a, b in
                        zip(tissue_tpm[t1], tissue_tpm[t2])]
            theset = math.sqrt(sum(theset))
            eucvalues.append(theset)
    eucvalues = np.array(eucvalues).reshape((c,c))
    return eucvalues

def manhattan(tissue_tpm):
    eucvalues = []
    for t1 in tissue_tpm:
        print ("Calculating %s set..." % t1)
        for t2 in tissue_tpm:
            theset = [(abs(a-b)) for a, b in
                        zip(tissue_tpm[t1], tissue_tpm[t2])]
            theset = sum(theset)
            eucvalues.append(theset)
    eucvalues = np.array(eucvalues).reshape((c,c))
    print (eucvalues)

manhattan(tissue_tpm)

# if __name__ == '__main__':
#     euclidean("../Data/HPA_processed.tsv")
