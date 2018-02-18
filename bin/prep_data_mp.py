#!/usr/bin/env python3
import math, sys, time, os
import numpy as np
import pandas as pd
# import multiprocessing as mp
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
c = len(genes)
reftype = str(sys.argv[2])

def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

eu_output_name = "../Data/{}_Euclidean_MP".format(reftype)
eu_output = open(eu_output_name+'.tsv', 'w')
eu_output.write('\t'.join(genes[:c])+'\n')
eu_mcl_name = "../Data/{}_EuclideanMCL_MP".format(reftype)
eu_mcl = open(eu_mcl_name+'.tsv', 'w')

def calc_matrix(enum, array):
    eu = []
    for j in range(c):
        vect2 = ar[j, ]
        eu.append(spd.euclidean(array, vect2))
    if (enum+1) % 25 == 0:
        print ("Normalising Euclidean {}/{}".format(enum+1, c))
    eu = normalise(eu)
    return [enum, genes[enum], eu]

results = Parallel(n_jobs=-1)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)

# def iterate(i):#, output):
#     eu = []
#     # vect1 = ar[i,]
#     for j in range(c):
#         # location = "[{}, {}]".format(i+1,j+1)
#         # print ("Position: {}".format(location))
#         vect2 = ar[j,]
#         eu.append(spd.euclidean(ar[i, ], vect2))
#     if (i+1) % 25 == 0:
#         print ("Normalising Euclidean {}/{}".format(i+1, c))
#     eu = normalise(eu)
#     return (i, genes[i], eu)
#     # output.put((i, genes[i], eu))
#
### Using python's multiprocessing library
# # output = mp.Queue()
# # processes = [mp.Process(target=iterate, args=(j, output)) for j in range(c)]
# # for p in processes:
# #     p.start()
# # for p in processes:
# #     p.join()
# # results = [output.get() for p in processes]
# # results.sort()
#
# # pool = mp.Pool(processes=8)
# # results = [pool.apply(iterate, args=(q, )) for q in range(c)]
# # print (results)
#
# results = Parallel(n_jobs=2)(delayed(iterate)(ar[i, ]) for i in range(c))
# for e in results:
#     print (e)
# # Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))]
### Which did not work out very well?

print ("Writing output...")
for r in results:
    eu_output.write(str(r[1])+'\t'+'\t'.join([str(k) for k in r[2]])+'\n')
    for h in range(int(r[0]), c):
        eu_mcl.write(str(r[1])+'\t'+genes[h]+'\t'+str(r[2][h])+'\n')

end = time.gmtime(time.time()-start)
print ("Distance matrix of {} for Euclidean metric completed in {}."
        .format(sys.argv[1], time.strftime("%Hh %Mm %Ss", end)))
