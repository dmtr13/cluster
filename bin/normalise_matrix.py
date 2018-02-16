#!/usr/bin/env python3
import sys, os, time
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, normalize

"""
This modular script takes in a generated distance matrix file and normalise it.
"""

print ("Processing matrix {}".format(sys.argv[1]))
start = time.time()

# pd.set_option("display.precision", 10)
df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)

def normalise(inlist):
    """
    Normalisation done by dividing each datapoint by the datapoint with the
    highest value in one particular set/list.
    """
    maxi = max(inlist)
    newlist = [i/maxi for i in inlist]
    return newlist

x_scaled = np.array([normalise(z) for z in df.values])
df = pd.DataFrame(x_scaled, columns=df.columns, index=df.index)

filename, file_ext = os.path.splitext(sys.argv[1])
outname = filename+"_Normalised"+file_ext
df.to_csv(path_or_buf=outname, sep='\t')

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
