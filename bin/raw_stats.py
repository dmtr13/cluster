#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, time, sys, math
import numpy as np

"""
Calculate for the input the number of genes that have zero expressions and not.
"""
filename_extention = os.path.basename(sys.argv[1])
directory = os.path.dirname(sys.argv[1])
fn, ext = os.path.splitext(filename_extention)

newfilename = directory+'/'+fn+".stats"
outstat = open(newfilename, 'w+')

count_genes, count_zeros = 0, 0
gene, gene_zero = [], []
with open(sys.argv[1], 'r') as proc:
    for z, line in enumerate(proc):
        if z == 0:
            line = line.split('\t')
            outstat.write("Total No. of Tissues:\t{}\n".format(len(line[1:])))
        else:
            count_genes += 1
            line = line.split()
            val = sum([float(i) for i in line[1:]])
            if val == 0:
                count_zeros += 1
                gene_zero.append(line[0])
            else:
                gene.append(line[0])
    outstat.write("Total No. of Genes:\t{}\n".format(str(count_genes)))
    outstat.write("Total No. of Genes with Non-zero Expression:\t{}\n".format(
                    str(count_genes-count_zeros)))
    outstat.write("Total No. of Genes with Zero Expression in All Tissues:\t{}\n".format(
                    str(count_zeros)))
    outstat.write("Genes with Non-zero Expression:\n")
    for g in gene:
        outstat.write(g+'\n')
    outstat.write("Genes with Zero Expression in All Tissues:\n")
    for gz in gene_zero:
        outstat.write(gz+'\n')

print ("Done.")
