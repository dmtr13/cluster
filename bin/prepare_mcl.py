#!/usr/bin/env python3
import sys, os
import pandas as pd

print ("Preparing {} for MCL...".format(sys.argv[1]))
"""
Takes a normalised matrix and converts into a 3-column format acceptable
for MCL.
"""

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
# Load the pandas DataFrame
# iterate through each row/line, appending into a list (or maybe a file, if
# it is more memory efficient) pairwise combination of the gene names
# and on the third column the values from the matrix.


filename, file_ext = os.path.splitext(sys.argv[1])
outname = filename+"_MCLin"+file_ext

print ("Done!")
