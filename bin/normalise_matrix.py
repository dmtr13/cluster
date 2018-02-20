#!/usr/bin/env python3
import math, sys, time, os, tempfile, shutil
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from joblib import Parallel, delayed

"""
This modular script takes in a generated distance matrix file and normalise it.
"""

print ("Processing matrix {}".format(sys.argv[1]))
start = time.time()

# pd.set_option("display.precision", 10)
df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
ar = np.array(df)
c = len(ar)

path = tempfile.mkdtemp()
outpath = os.path.join(path, 'matrix_temp.mmap')
outfile = np.memmap(outpath, dtype=float, shape=df.shape, mode='w+')

def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

def norm_matrix(enum, array):
    if (enum+1) % 100 == 0:
        print ("CHECKPOINT {}/{}...".format(enum+1, c))
    else:
        print ("Processing {}/{}...".format(enum+1, c))
    outfile[enum, ] = normalise(array)
    return True

results = Parallel(n_jobs=-1)(delayed(norm_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)

# x_scaled = np.array([normalise(z) for z in df.values])
df_N = pd.DataFrame(outfile, columns=df.columns, index=df.index)

filename, file_ext = os.path.splitext(sys.argv[1])
outname = filename+"_Normalised"+file_ext
df_N.to_csv(path_or_buf=outname, sep='\t')

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
