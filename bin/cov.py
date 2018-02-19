#!/usr/bin/env python3
import math, sys, time, os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a normalised matrix for Euclidean and Manhattan metrics and also
the three-column format accepted by MCL. Multiprocessing enabled.
"""
start = time.time()
print ("Creating distance matrices using Euclidean metric...")

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = 50#len(genes)
reftype = str(sys.argv[2])


def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

# eu_output_name = "../Data/{}_Euclidean_MP".format(reftype)
# eu_output = open(eu_output_name+'.tsv', 'w')
# eu_output.write('\t'.join(genes[:c])+'\n')
# eu_mcl_name = "../Data/{}_EuclideanMCL_MP".format(reftype)
# eu_mcl = open(eu_mcl_name+'.tsv', 'w')

covar = np.empty([c,c], dtype=float)

def calc_matrix(enum, array): ## COVARIANCE
    cv = []
    if (enum + 1) % 2 == 0:
        print ("Processing {}/{}...".format(enum+1, c))
    for j in range(c):
        vect2 = ar[j, ]
        covariance = np.cov(array, vect2, bias=True)
        # if enum == 0 and j == 0:
        print (enum, j, covariance[0,0])
        covar[enum, j] == covariance[0,0]
    #     cv.append(covariance[1,0])
    # covar[enum, ] = cv
    return True
#
results = Parallel(n_jobs=-1)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)
#
# print ("Writing output...")
# for r in results:
#     eu_output.write(str(r[1])+'\t'+'\t'.join([str(k) for k in r[2]])+'\n')
#     for h in range(int(r[0]), c):
#         eu_mcl.write(str(r[1])+'\t'+genes[h]+'\t'+str(r[2][h])+'\n')
#
# end = time.gmtime(time.time()-start)
# print ("Distance matrix of {} for Euclidean metric completed in {}."
#         .format(sys.argv[1], time.strftime("%Hh %Mm %Ss", end)))

## https://stackoverflow.com/questions/45487145/pandas-correlation-between-list-of-columns-x-whole-dataframe
