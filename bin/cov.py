#!/usr/bin/env python3
import math, sys, time, os, tempfile, shutil
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a normalised matrix for Euclidean and Manhattan metrics and also
the three-column format accepted by MCL. Multiprocessing enabled.
"""
start = time.time()
print ("Creating a partial correlation via inverse covariance matrix...")

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = 50#len(genes)

if "HPA" in str(sys.argv[1]):
    reftype = "HPA"
elif "GTEx" in str(sys.argv[1]):
    reftype = "GTEx"


def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

eu_output_name = "../Data/{}_InCov".format(reftype)
eu_output = open(eu_output_name+'.tsv', 'w')
# eu_output.write('\t'.join(genes[:c])+'\n')
eu_mcl_name = "../Data/{}_EuclideanMCL_MP".format(reftype)
# eu_mcl = open(eu_mcl_name+'.tsv', 'w')

path = tempfile.mkdtemp()
covar_path = os.path.join(path, 'incov.mmap')
covar = np.memmap(covar_path, dtype=float, shape=(c,c), mode='w+')

def calc_matrix(enum, array): ## COVARIANCE
    if (enum+1) % 100 == 0:
        print ("CHECKPOINT {}/{}...".format(enum+1, c))
    else:
        print ("Processing {}/{}...".format(enum+1, c))
    temp = []
    for j in range(c):
        vect2 = ar[j, ]
        covariance = np.cov(array, vect2, bias=True)
        # print (enum, j, covariance)
        # covar[enum, j] == covariance[1,0]
        temp.append(covariance[1,0])
    covar[enum, ] = temp
    return True
#
results = Parallel(n_jobs=-2)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)

print ("Inversing covariance matrix...")
incov = np.linalg.inv(covar)
incov = np.array([normalise(x) for x in incov])
try:
    shutil.rmtree(path)
    print ("Temporary file deleted.")
except:
    pass

df_incov = pd.DataFrame(incov, columns=genes[:c], index=genes[:c])
df_incov.to_csv(path_or_buf=eu_output_name+'.tsv', sep='\t')
print (df_incov.iloc[:3, :3])
# with open(eu_mcl_name+'.tsv', 'w') as o:
#     for z in range(c):


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
## https://stackoverflow.com/questions/34140560/accessing-and-altering-a-global-array-using-python-joblib
