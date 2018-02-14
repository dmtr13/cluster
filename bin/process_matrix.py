#!/usr/bin/env python3
import sys
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
# arr = np.loadtxt(strip_first_col(sys.argv[1]), skiprows=1)
# arr = normalize(arr, axis=1, norm='l1')
# print (arr)

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
# print (df.iloc[[0]])
for i in range(len(df)):
    maxval = (df.iloc[[i]]).idmax()
    print (maxval)
    # print (df.iloc[[i]])


# scaler = MinMaxScaler()
# scaled = scaler.fit_transform(df[1])
# print (scaled[0])
# df.loc[:,:] = scaled
# print (df[0])
