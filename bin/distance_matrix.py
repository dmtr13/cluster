#!/bin/usr/env python3

def euclidean(filename):
    print ("a")
    filename = open(filename, 'r')
    filename = filename.read().splitlines()
    header = filename[0].split()
    dimension = len(header)
    genes = []
    values = []
    indexer = 0
    for lines in filename[1:]:
        lines = lines.split()
        genes.append(lines[0])
        values.append(lines[1:])
    print (header)

if __name__ == '__main__':
    euclidean("../Data/HPA_processed.tsv")
