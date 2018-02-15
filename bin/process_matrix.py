#!/usr/bin/env python3
import sys, os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, normalize

print ("Processing matrix {}".format(sys.argv[1]))

def strip_first_col(fname, delimiter='\t'):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue
# pd.set_option("display.precision", 10)
df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)

def normalise(inlist):
    maxi = max(inlist)
    newlist = [i/maxi for i in inlist]
    return newlist
x = df.values
x_scaled = np.array([normalise(z) for z in x])
df = pd.DataFrame(x_scaled, columns=df.columns, index=df.index)

filename, file_ext = os.path.splitext(sys.argv[1])
outname = filename+"_Normalised"+file_ext
df.to_csv(path_or_buf=outname, sep='\t')

print ("Done!")
