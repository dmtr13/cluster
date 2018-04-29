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
top10 = math.ceil(c*0.9) ## 90% cut-off threshold
reftype = str(sys.argv[2])

def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

def rel_euclidean(vect1, vect2):
    summa = 0
    for v1, v2 in zip(vect1, vect2):
        ## The error is because it's dividing by 0. E.g. v1=0, v2=0, which
        ## then makes it nan.
        distance = (v1-v2)**2
        divider = (math.fabs(v1) + math.fabs(v2))/2
        if divider == 0:
            summa += 0
        else:
            summa += distance/divider
    return math.sqrt(summa)

def prune_top10pc(vect):
    """ Takes the values of the top 10%, everything else will be 0.
    """
    from scipy.stats import rankdata as rd
    ranking = [int(x-1) for x in rd(vect)]
    ranked_vect = []
    for x in ranking:
        if x >= top10:
            ranked_vect.extend([vect[x]])
        else:
            ranked_vect.extend([0])
    return ranked_vect

def normalise_by_max(inlist):
    return [i-max(inlist) for i in inlist]

eu_output_name = "../Data/{}_RelEuclidean".format(reftype)
eu_output = open(eu_output_name+'.tsv', 'w')
eu_output.write('\t'.join(genes[:c])+'\n')
eu_mcl_name = "../Data/{}_RelEuclideanMCL".format(reftype)
eu_mcl = open(eu_mcl_name+'.tsv', 'w')

def calc_matrix(enum, array):
    # print ("Processing {}/{}...".format(enum, c))
    eu = []
    for j in range(c):
        vect2 = ar[j, ]
        eu.append(rel_euclidean(array, vect2))
    # print (enum, eu)
    eu = [max(eu) - el for el in eu]
    eu = prune_top10pc(eu)
    # print (enum, eu)
    if (enum+1) % 25 == 0:
        print ("Checkpoint: {}/{}".format(enum+1, c))
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
